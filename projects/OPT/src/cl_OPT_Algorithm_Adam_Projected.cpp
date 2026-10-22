/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Algorithm_Adam_Projected.cpp
 *
 * Adam with an exact one-constraint projection (see header). MORIS owns the loop, as in Algorithm_MMA.
 */

#include <algorithm>
#include <cmath>
#include <sstream>

#include "cl_OPT_Algorithm_Adam_Projected.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

namespace moris::opt
{
    //------------------------------------------------------------------------------

    std::vector< real >
    Algorithm_Adam_Projected::parse_pattern( const std::string& aString )
    {
        std::vector< real > tOut;
        std::stringstream   tStream( aString );
        std::string         tToken;
        while ( std::getline( tStream, tToken, ',' ) )
        {
            if ( !tToken.empty() )
            {
                tOut.push_back( std::stod( tToken ) );
            }
        }
        return tOut;
    }

    //------------------------------------------------------------------------------

    Algorithm_Adam_Projected::Algorithm_Adam_Projected( const Parameter_List& aParameterList )
            : mLearningRate( aParameterList.get< real >( "learning_rate" ) )
            , mRatePattern( parse_pattern( aParameterList.get< std::string >( "learning_rate_pattern" ) ) )
            , mDecayFinal( aParameterList.get< real >( "lr_decay_final" ) )
            , mBeta1( aParameterList.get< real >( "adam_beta1" ) )
            , mBeta2( aParameterList.get< real >( "adam_beta2" ) )
            , mEpsilon( aParameterList.get< real >( "adam_epsilon" ) )
            , mMaskPattern( parse_pattern( aParameterList.get< std::string >( "projection_mask" ) ) )
            , mDirection( aParameterList.get< std::string >( "projection_direction" ) )
            , mMode( aParameterList.get< std::string >( "projection_mode" ) )
            , mProjectionTol( aParameterList.get< real >( "projection_tol" ) )
            , mProjectionMaxIts( aParameterList.get< sint >( "projection_max_its" ) )
            , mTangentMoments( aParameterList.get< bool >( "tangent_moments" ) )
    {
        mRestartIndex         = aParameterList.get< sint >( "restart_index" );
        mMaxIterationsInitial = aParameterList.get< sint >( "max_its" );
        mMaxIterations        = mMaxIterationsInitial;

        MORIS_ERROR( mDirection == "uniform" || mDirection == "gradient",
                "Algorithm_Adam_Projected: projection_direction must be \"uniform\" or \"gradient\"" );
        MORIS_ERROR( mMode == "equality" || mMode == "inequality",
                "Algorithm_Adam_Projected: projection_mode must be \"equality\" or \"inequality\"" );
        MORIS_ERROR( mDecayFinal > 0.0, "Algorithm_Adam_Projected: lr_decay_final must be > 0" );

        MORIS_LOG_INFO( "Adam-projected initialized: max_its=%d, learning_rate=%.4e, lr_decay_final=%.3e, direction=%s, mode=%s",
                (int)mMaxIterations, mLearningRate, mDecayFinal, mDirection.c_str(), mMode.c_str() );
    }

    //------------------------------------------------------------------------------

    Algorithm_Adam_Projected::~Algorithm_Adam_Projected() {}

    //------------------------------------------------------------------------------

    uint
    Algorithm_Adam_Projected::solve( uint aCurrentOptAlgInd, std::shared_ptr< Problem > aOptProb )
    {
        Tracer tTracer( "OptimizationAlgorithm", "AdamP", "Solve" );

        mRunning          = opt::Task::wait;
        mCurrentOptAlgInd = aCurrentOptAlgInd;
        mProblem          = aOptProb;
        mNumForwardSolves = 0;

        if ( mRestartIndex > 0 )
        {
            gLogger.set_opt_iteration( mRestartIndex );
        }

        if ( par_rank() == 0 )
        {
            this->adam_solve();

            mRunning = opt::Task::exit;
            this->communicate_running_status();
        }
        else
        {
            this->dummy_solve();
        }

        uint tOptIter = gLogger.get_opt_iteration();
        gLogger.set_iteration( "OPT", "Manager", "Perform", tOptIter );

        return tOptIter;
    }

    //------------------------------------------------------------------------------

    void
    Algorithm_Adam_Projected::evaluate( const std::vector< real >& aX, real& aF, real& aG )
    {
        Vector< real > tADVs( aX.size() );
        for ( uint i = 0; i < aX.size(); ++i )
        {
            tADVs( i ) = aX[ i ];
        }
        this->compute_design_criteria( tADVs );
        ++mNumForwardSolves;
        aF = this->get_objectives()( 0 );
        aG = mProblem->get_num_constraints() > 0 ? this->get_constraints()( 0 ) : 0.0;
    }

    //------------------------------------------------------------------------------

    void
    Algorithm_Adam_Projected::adam_solve()
    {
        const uint tN    = mProblem->get_num_advs();
        const uint tNCon = mProblem->get_num_constraints();
        MORIS_ERROR( tNCon <= 1, "Algorithm_Adam_Projected supports at most one constraint (got %d)", (int)tNCon );

        std::vector< real > tX( tN ), tLo( tN ), tHi( tN ), tRate( tN, 1.0 ), tMask( tN, 1.0 );
        for ( uint i = 0; i < tN; ++i )
        {
            tLo[ i ] = mProblem->get_lower_bounds()( i );
            tHi[ i ] = mProblem->get_upper_bounds()( i );
            tX[ i ]  = std::min( std::max( mProblem->get_advs()( i ), tLo[ i ] ), tHi[ i ] );
            if ( !mRatePattern.empty() ) tRate[ i ] = mRatePattern[ i % mRatePattern.size() ];
            if ( !mMaskPattern.empty() ) tMask[ i ] = mMaskPattern[ i % mMaskPattern.size() ];
        }
        auto tClip = [ & ]( std::vector< real >& aV ) {
            for ( uint i = 0; i < tN; ++i ) aV[ i ] = std::min( std::max( aV[ i ], tLo[ i ] ), tHi[ i ] );
        };

        std::vector< real > tM( tN, 0.0 ), tV( tN, 0.0 ), tGf( tN ), tGc( tN, 0.0 ), tD( tN ), tTrial( tN );

        real tF = 0.0, tG = 0.0;
        this->evaluate( tX, tF, tG );

        const real tDecayRate = ( mMaxIterations > 1 ) ? std::log( mDecayFinal ) / (real)( mMaxIterations - 1 ) : 0.0;

        for ( uint tIter = 0; tIter < (uint)mMaxIterations; ++tIter )
        {
            // gradients at the current (projected) iterate; the iterate is accepted
            Vector< real > tADVs( tN );
            for ( uint i = 0; i < tN; ++i ) tADVs( i ) = tX[ i ];
            this->compute_design_criteria_gradients( tADVs );
            const Matrix< DDRMat >& tDf = this->get_objective_gradients();
            for ( uint i = 0; i < tN; ++i ) tGf[ i ] = tDf( 0, i );
            if ( tNCon == 1 )
            {
                const Matrix< DDRMat >& tDg = this->get_constraint_gradients();
                for ( uint i = 0; i < tN; ++i ) tGc[ i ] = tDg( 0, i );
            }
            this->write_advs_to_file( tADVs );

            // projection direction d (masked) and constraint slope gc . d
            for ( uint i = 0; i < tN; ++i ) tD[ i ] = tMask[ i ] * ( mDirection == "gradient" ? tGc[ i ] : 1.0 );
            real tSlope = 0.0, tDD = 0.0, tFD = 0.0;
            for ( uint i = 0; i < tN; ++i )
            {
                tSlope += tGc[ i ] * tD[ i ];
                tDD += tD[ i ] * tD[ i ];
            }
            const bool tProject = ( tNCon == 1 ) && std::abs( tSlope ) > 0.0;

            // Adam moments on the objective gradient, optionally with its component along the constraint normal
            // removed (the projection undoes that component anyway; keeping it would make the momentum fight it)
            std::vector< real > tGt = tGf;
            if ( mTangentMoments && tProject && ( mMode == "equality" || tG > -mProjectionTol ) )
            {
                real tGG = 0.0;
                for ( uint i = 0; i < tN; ++i )
                {
                    tFD += tGf[ i ] * tGc[ i ];
                    tGG += tGc[ i ] * tGc[ i ];
                }
                if ( tGG > 0.0 )
                    for ( uint i = 0; i < tN; ++i ) tGt[ i ] -= tFD / tGG * tGc[ i ];
            }

            const real tLR = mLearningRate * std::exp( tDecayRate * (real)tIter );
            const real tB1 = 1.0 - std::pow( mBeta1, (real)( tIter + 1 ) );
            const real tB2 = 1.0 - std::pow( mBeta2, (real)( tIter + 1 ) );
            real       tStepMax = 0.0;
            for ( uint i = 0; i < tN; ++i )
            {
                tM[ i ]     = mBeta1 * tM[ i ] + ( 1.0 - mBeta1 ) * tGt[ i ];
                tV[ i ]     = mBeta2 * tV[ i ] + ( 1.0 - mBeta2 ) * tGt[ i ] * tGt[ i ];
                real tStep  = -tLR * tRate[ i ] * ( tM[ i ] / tB1 ) / ( std::sqrt( tV[ i ] / tB2 ) + mEpsilon );
                tTrial[ i ] = tX[ i ] + tStep;
                tStepMax    = std::max( tStepMax, std::abs( tStep ) );
            }
            tClip( tTrial );

            // projection: x(t) = trial + t d; linearized predictor g(x(t)) ~ g_k + gc.(trial - x_k) + t gc.d
            real tT = 0.0;
            if ( tProject )
            {
                real tPred = tG;
                for ( uint i = 0; i < tN; ++i ) tPred += tGc[ i ] * ( tTrial[ i ] - tX[ i ] );
                if ( mMode == "equality" || tPred > 0.0 ) tT = -tPred / tSlope;
            }
            std::vector< real > tXNew( tN );
            auto                tAt = [ & ]( real aT ) {
                for ( uint i = 0; i < tN; ++i ) tXNew[ i ] = tTrial[ i ] + aT * tD[ i ];
                tClip( tXNew );
            };
            tAt( tT );
            real tFNew = 0.0, tGNew = 0.0;
            this->evaluate( tXNew, tFNew, tGNew );

            // secant corrections on the true constraint
            int  tCorr  = 0;
            real tTPrev = tT, tGPrev = tGNew, tSec = tSlope;
            auto tNeeds = [ & ]( real aG ) {
                return tProject && ( mMode == "equality" ? std::abs( aG ) > mProjectionTol : aG > mProjectionTol );
            };
            while ( tNeeds( tGNew ) && tCorr < mProjectionMaxIts )
            {
                real tTNext = tTPrev - tGPrev / tSec;
                tAt( tTNext );
                this->evaluate( tXNew, tFNew, tGNew );
                ++tCorr;
                if ( std::abs( tTNext - tTPrev ) > 0.0 )
                {
                    real tS = ( tGNew - tGPrev ) / ( tTNext - tTPrev );
                    if ( std::abs( tS ) > 1.0e-3 * std::abs( tSlope ) ) tSec = tS;    // keep a sane slope
                }
                tTPrev = tTNext;
                tGPrev = tGNew;
            }

            tX = tXNew;
            tF = tFNew;
            tG = tGNew;

            MORIS_LOG_INFO( "ADAMP it=%d fsolves=%d obj=%.8e con=%.4e lr=%.4e step_max=%.4e t=%.4e corrections=%d",
                    (int)tIter + 1, (int)mNumForwardSolves, tF, tG, tLR, tStepMax, tTPrev, tCorr );
        }

        // final gradients so the Problem is consistent at the last design (and its restart file is written)
        Vector< real > tADVs( tN );
        for ( uint i = 0; i < tN; ++i ) tADVs( i ) = tX[ i ];
        this->write_advs_to_file( tADVs );
    }

    //------------------------------------------------------------------------------

}    // namespace moris::opt

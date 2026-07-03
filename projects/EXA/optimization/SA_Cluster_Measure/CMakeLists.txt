#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#
moris_add_example(
    NAME       SA_Cluster_Measure
    TEST_BASE  SA_Cluster_Measure    # ctest name stays SA_Cluster_Measure-1-procs
    PROCS      1
    INPUTS     Cluster_Measure_2Mat_SA    # deck basename differs from TEST_BASE
)

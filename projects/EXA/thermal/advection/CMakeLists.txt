#
# Copyright (c) 2022 University of Colorado
# Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
#
#------------------------------------------------------------------------------------
#

add_subdirectory(Channel_with_Four_Cylinders_Static)

# Disconnected: these decks use the pre-modern ParameterList API and no longer compile.
# The leaves are already migrated to moris_add_example(); reconnect by uncommenting
# once the decks are ported to Module_Parameter_Lists (see EXA_RUNNER_RFC.md).
# add_subdirectory(Channel_with_Four_Cylinders_Transient)
# add_subdirectory(Two_Channels_with_Separation_Wall_Transient)


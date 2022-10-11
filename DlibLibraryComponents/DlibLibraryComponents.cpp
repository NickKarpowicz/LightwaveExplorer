// DlibLibraryComponents.cpp : Defines the functions for the static library.
//

#include "pch.h"
#include "framework.h"

// Build the necessary components of dlib for optimization
// This is part of the Dlib library, whose repo should 
// be cloned to ../dlib
// See dlib.net for more info
// Copyright (C) 2006  Davis E. King (davis@dlib.net)
// License: Boost Software License   See LICENSE.txt for the full license.
#include "dlib/test_for_odr_violations.cpp"
#include "dlib/threads/thread_pool_extension.cpp"
#include "dlib/global_optimization/global_function_search.cpp"
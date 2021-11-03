/usr/local/apps/matlab/bin/mex CXXFLAGS="\$CXXFLAGS -std=c++11 -fopenmp" CXXOPTIMFLAGS='\$CXXOPTIMFLAGS -Ofast -DNDEBUG' LDOPTIMFLAGS="$LDOPTIMFLAGS -fopenmp -O2" -I"." fwCoreMex.cpp
/usr/local/apps/matlab/bin/mex CXXFLAGS="\$CXXFLAGS -std=c++11 -fopenmp" CXXOPTIMFLAGS='\$CXXOPTIMFLAGS -Ofast -DNDEBUG' LDOPTIMFLAGS="$LDOPTIMFLAGS -fopenmp -O2" -I"." fwCoreSparseMex.cpp


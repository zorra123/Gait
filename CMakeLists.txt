cmake_minimum_required (VERSION 2.6)



project (walk)

include_directories (include/)


set (SRC 
        #include/MotionModule.h
        include/vector.h
		include/matrix.h
        #include/MX28.h
        #include/MotionStatus.h
		#include/CM730.h
        #include/FSR.h
        #include/JointData.h
        #include/MotionManager.h
		include/point.h
        include/walking.h
        src/vector.cpp
		src/matrix.cpp
        #src/MX28.cpp
        #src/MotionStatus.cpp
		#src/CM730.cpp
        #src/JointData.cpp
        #src/MotionManager.cpp
		src/point.cpp
		include/defines.h
        src/walking.cpp)

add_library (walking ${SRC} )

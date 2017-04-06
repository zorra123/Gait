cmake_minimum_required (VERSION 2.6)



project (walk)

include_directories (include/)


set (SRC 
        #include/MotionModule.h
	include/Vector.h 
	include/Matrix.h 
        #include/MX28.h
        #include/MotionStatus.h
	include/Kinematics.h
        #include/CM730.h
        #include/FSR.h
        #include/JointData.h
        #include/MotionManager.h
	include/Point.h
	Walking.h  
	src/Vector.cpp 
	src/Matrix.cpp 
        #src/MX28.cpp
        #src/MotionStatus.cpp
	src/Kinematics.cpp
        #src/CM730.cpp
        #src/JointData.cpp
        #src/MotionManager.cpp
	src/Point.cpp
        include/const.h
	Walking.cpp)

add_library (walking ${SRC} )

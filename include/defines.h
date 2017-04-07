
namespace darwin {
    const int CENTER_VALUE = 2048;
    const double TIME_UNIT = 8;

    namespace joint {
        enum joint_id {
            ID_R_SHOULDER_PITCH = 1,
            ID_L_SHOULDER_PITCH = 2,
            ID_R_SHOULDER_ROLL = 3,
            ID_L_SHOULDER_ROLL = 4,
            ID_R_ELBOW = 5,
            ID_L_ELBOW = 6,
            ID_R_HIP_YAW = 7,
            ID_L_HIP_YAW = 8,
            ID_R_HIP_ROLL = 9,
            ID_L_HIP_ROLL = 10,
            ID_R_HIP_PITCH = 11,
            ID_L_HIP_PITCH = 12,
            ID_R_KNEE = 13,
            ID_L_KNEE = 14,
            ID_R_ANKLE_PITCH = 15,
            ID_L_ANKLE_PITCH = 16,
            ID_R_ANKLE_ROLL = 17,
            ID_L_ANKLE_ROLL = 18,
            ID_HEAD_PAN = 19,
            ID_HEAD_TILT = 20,
            NUMBER_OF_JOINTS
        };

        const double CW_LIMIT_R_SHOULDER_ROLL = -75.0; // degree
        const double CCW_LIMIT_R_SHOULDER_ROLL = 135.0; // degree
        const double CW_LIMIT_L_SHOULDER_ROLL = -135.0; // degree
        const double CCW_LIMIT_L_SHOULDER_ROLL = 75.0; // degree
        const double CW_LIMIT_R_ELBOW = -95.0; // degree
        const double CCW_LIMIT_R_ELBOW = 70.0; // degree
        const double CW_LIMIT_L_ELBOW = -70.0; // degree
        const double CCW_LIMIT_L_ELBOW = 95.0; // degree
        const double CW_LIMIT_R_HIP_YAW = -123.0; // degree
        const double CCW_LIMIT_R_HIP_YAW = 53.0; // degree
        const double CW_LIMIT_L_HIP_YAW = -53.0; // degree
        const double CCW_LIMIT_L_HIP_YAW = 123.0; // degree
        const double CW_LIMIT_R_HIP_ROLL = -45.0; // degree
        const double CCW_LIMIT_R_HIP_ROLL = 59.0; // degree
        const double CW_LIMIT_L_HIP_ROLL = -59.0; // degree
        const double CCW_LIMIT_L_HIP_ROLL = 45.0; // degree
        const double CW_LIMIT_R_HIP_PITCH = -100.0; // degree
        const double CCW_LIMIT_R_HIP_PITCH = 29.0; // degree
        const double CW_LIMIT_L_HIP_PITCH = -29.0; // degree
        const double CCW_LIMIT_L_HIP_PITCH = 100.0; // degree
        const double CW_LIMIT_R_KNEE = -6.0; // degree
        const double CCW_LIMIT_R_KNEE = 130.0; // degree
        const double CW_LIMIT_L_KNEE = -130.0; // degree
        const double CCW_LIMIT_L_KNEE = 6.0; // degree
        const double CW_LIMIT_R_ANKLE_PITCH = -72.0; // degree
        const double CCW_LIMIT_R_ANKLE_PITCH = 80.0; // degree
        const double CW_LIMIT_L_ANKLE_PITCH = -80.0; // degree
        const double CCW_LIMIT_L_ANKLE_PITCH = 72.0; // degree
        const double CW_LIMIT_R_ANKLE_ROLL = -44.0; // degree
        const double CCW_LIMIT_R_ANKLE_ROLL = 63.0; // degree
        const double CW_LIMIT_L_ANKLE_ROLL = -63.0; // degree
        const double CCW_LIMIT_L_ANKLE_ROLL = 44.0; // degree
        const double CW_LIMIT_HEAD_PAN = -90.0; // degree
        const double CCW_LIMIT_HEAD_PAN = 90.0; // degree
        const double CW_LIMIT_HEAD_TILT = -25.0; // degree
        const double CCW_LIMIT_HEAD_TILT = 55.0; // degree
    }

    namespace kinematics {
        const double CAMERA_DISTANCE = 33.2; //mm
        const double EYE_TILT_OFFSET_ANGLE = 40.0; //degree
        const double LEG_SIDE_OFFSET = 37.0; //mm
        const double THIGH_LENGTH = 93.0; //mm
        const double CALF_LENGTH = 93.0; //mm
        const double ANKLE_LENGTH = 33.5; //mm
        const double LEG_LENGTH = 219.5; //mm (THIGH_LENGTH + CALF_LENGTH + ANKLE_LENGTH)
    };
}

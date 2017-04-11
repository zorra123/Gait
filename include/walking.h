/*
 *   Walking.h
 *
 *   Author: ROBOTIS
 *
 */

#ifndef _WALKING_ENGINE_H_
#define _WALKING_ENGINE_H_

#include <string.h>

namespace darwin {
    class walking {
    public:
        walking();

        // walking initial pose
        double x_offset;
        double y_offset;
        double z_offset;
        double a_offset;
        double p_offset;
        double r_offset;

        // walking control
        double period_time;
        double dsp_ratio;
        double step_fb_ratio;
        double x_move_amplitude;
        double y_move_amplitude;
        double z_move_amplitude;
        double a_move_amplitude;
        bool a_move_aim_on;

        // Balance control
        bool m_balance_enable;
        double balance_knee_gain;
        double balance_ankle_pitch_gain;
        double balance_hip_roll_gain;
        double balance_ankle_roll_gain;
        double y_swap_amplitude;
        double z_swap_amplitude;
        double arm_swing_gain;
        double pelvis_offset;
        double hip_pitch_offset;

        int p_gain;
        int i_gain;
        int d_gain;

        ~walking();

        void initialize();

        void start();

        void stop();

        void process();

        bool is_running();

        void set_p_gain(int id, int pgain);

        void set_angle(int id, double angle);
/*		todo
		void LoadINISettings(minIni* ini);
		void LoadINISettings(minIni* ini, const std::string &section);
		void SaveINISettings(minIni* ini);
		void SaveINISettings(minIni* ini, const std::string &section);
	*/
    private:
        enum local_joint_id {
            R_HIP_YAW,
            R_HIP_ROLL,
            R_HIP_PITCH,
            R_KNEE,
            R_ANKLE_PITCH,
            R_ANKLE_ROLL,
            L_HIP_YAW,
            L_HIP_ROLL,
            L_HIP_PITCH,
            L_KNEE,
            L_ANKLE_PITCH,
            L_ANKLE_ROLL,
            R_ARM_SWING,
            L_ARM_SWING,
            NUM_OF_JOINTS
        };

        enum position {
            R_X,
            R_Y,
            R_Z,
            R_ROLL,
            R_PITCH,
            R_YAW,
            L_X,
            L_Y,
            L_Z,
            L_ROLL,
            L_PITCH,
            L_YAW,
        };

        enum phase {
            PHASE0 = 0,
            PHASE1 = 1,
            PHASE2 = 2,
            PHASE3 = 3
        };

        double m_period_time;
        double m_dsp_ratio;
        double m_ssp_ratio;
        double m_x_swap_period_time;
        double m_x_move_period_time;
        double m_y_swap_period_time;
        double m_y_move_period_time;
        double m_z_swap_period_time;
        double m_z_move_period_time;
        double m_a_move_period_time;
        double m_ssp_time;
        double m_ssp_time_start_l;
        double m_ssp_time_end_l;
        double m_ssp_time_start_r;
        double m_ssp_time_end_r;
        double m_phase_time1;
        double m_phase_time2;
        double m_phase_time3;

        double m_x_offset;
        double m_y_offset;
        double m_z_offset;
        double m_r_offset;
        double m_p_offset;
        double m_a_offset;

        double m_x_swap_phase_shift;
        double m_x_swap_amplitude;
        double m_x_swap_amplitude_shift;
        double m_x_move_phase_shift;
        double m_x_move_amplitude;
        double m_x_move_amplitude_shift;
        double m_y_swap_phase_shift;
        double m_y_swap_amplitude;
        double m_y_swap_amplitude_shift;
        double m_y_move_phase_shift;
        double m_y_move_amplitude;
        double m_y_move_amplitude_shift;
        double m_z_swap_phase_shift;
        double m_z_swap_amplitude;
        double m_z_swap_amplitude_shift;
        double m_z_move_phase_shift;
        double m_z_move_amplitude;
        double m_z_move_amplitude_shift;
        double m_a_move_phase_shift;
        double m_a_move_amplitude;
        double m_a_move_amplitude_shift;

        double m_pelvis_offset;
        double m_pelvis_swing;
        double m_hip_pitch_offset;
        double m_arm_swing_gain;

        bool m_ctrl_running;
        bool m_real_running;
        double m_time;

        double m_body_swing_y;
        double m_body_swing_z;

        static constexpr double wsin(double time, double period, double period_shift, double mag, double mag_shift);

        static bool compute_ik(double* out, double x, double y, double z, double a, double b, double c);

        void update_param_time();

        void update_param_move();

        void update_param_balance();
    };
}

#endif

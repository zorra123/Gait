/*
 *   Walking.cpp
 *
 *   Author: ROBOTIS
 *
 */

#include <iostream>
#include <math.h>
#include "vector.h"
#include "matrix.h"
#include "walking.h"
#include "defines.h"

using namespace darwin;


#define PI (3.14159265)


walking::walking() {
    x_offset = -10;
    y_offset = 5;
    z_offset = 20;
    r_offset = 0;
    p_offset = 0;
    a_offset = 0;
    hip_pitch_offset = 13.0;
    period_time = 600;
    dsp_ratio = 0.1;
    step_fb_ratio = 0.28;
    z_move_amplitude = 40;
    y_swap_amplitude = 20.0;
    z_swap_amplitude = 5;
    pelvis_offset = 3.0;
    arm_swing_gain = 1.5;
    balance_knee_gain = 0.3;
    balance_ankle_pitch_gain = 0.9;
    balance_hip_roll_gain = 0.5;
    balance_ankle_roll_gain = 1.0;
/* todo GAIN
    p_gain = JointData::P_GAIN_DEFAULT;//todo
    i_gain = JointData::I_GAIN_DEFAULT;//todo
    d_gain = JointData::D_GAIN_DEFAULT;//todo
*/

    x_move_amplitude = 0;
    y_move_amplitude = 0;
    a_move_amplitude = 0;
    a_move_aim_on = false;
    m_balance_enable = true;
/* todo set_angle */
    set_angle(joint::ID_R_SHOULDER_PITCH, -48.345);
    set_angle(joint::ID_L_SHOULDER_PITCH, 41.313);
    set_angle(joint::ID_R_SHOULDER_ROLL, -17.873);
    set_angle(joint::ID_L_SHOULDER_ROLL, 17.580);
    set_angle(joint::ID_R_ELBOW, 29.300);
    set_angle(joint::ID_L_ELBOW, -29.593);

    set_angle(joint::ID_HEAD_TILT, kinematics::EYE_TILT_OFFSET_ANGLE);

/* todo SetPGain*/
    set_p_gain(joint::ID_R_SHOULDER_PITCH, 8);
    set_p_gain(joint::ID_L_SHOULDER_PITCH, 8);
    set_p_gain(joint::ID_R_SHOULDER_ROLL, 8);
    set_p_gain(joint::ID_L_SHOULDER_ROLL, 8);
    set_p_gain(joint::ID_R_ELBOW, 8);
    set_p_gain(joint::ID_L_ELBOW, 8);

}


walking::~walking() {
}


/* todo INI
void walking::LoadINISettings(minIni* ini)
{
    LoadINISettings(ini, WALKING_SECTION);
}
void walking::LoadINISettings(minIni* ini, const std::string &section)
{
    double value = INVALID_VALUE;
    if((value = ini->getd(section, "x_offset", INVALID_VALUE)) != INVALID_VALUE)                x_offset = value;
    if((value = ini->getd(section, "y_offset", INVALID_VALUE)) != INVALID_VALUE)                y_offset = value;
    if((value = ini->getd(section, "z_offset", INVALID_VALUE)) != INVALID_VALUE)                z_offset = value;
    if((value = ini->getd(section, "roll_offset", INVALID_VALUE)) != INVALID_VALUE)             r_offset = value;
    if((value = ini->getd(section, "pitch_offset", INVALID_VALUE)) != INVALID_VALUE)            p_offset = value;
    if((value = ini->getd(section, "yaw_offset", INVALID_VALUE)) != INVALID_VALUE)              a_offset = value;
    if((value = ini->getd(section, "hip_pitch_offset", INVALID_VALUE)) != INVALID_VALUE)        hip_pitch_offset = value;
    if((value = ini->getd(section, "period_time", INVALID_VALUE)) != INVALID_VALUE)             period_time = value;
    if((value = ini->getd(section, "dsp_ratio", INVALID_VALUE)) != INVALID_VALUE)               dsp_ratio = value;
    if((value = ini->getd(section, "step_forward_back_ratio", INVALID_VALUE)) != INVALID_VALUE) step_fb_ratio = value;
    if((value = ini->getd(section, "foot_height", INVALID_VALUE)) != INVALID_VALUE)             z_move_amplitude = value;
    if((value = ini->getd(section, "swing_right_left", INVALID_VALUE)) != INVALID_VALUE)        y_swap_amplitude = value;
    if((value = ini->getd(section, "swing_top_down", INVALID_VALUE)) != INVALID_VALUE)          z_swap_amplitude = value;
    if((value = ini->getd(section, "pelvis_offset", INVALID_VALUE)) != INVALID_VALUE)           pelvis_offset = value;
    if((value = ini->getd(section, "arm_swing_gain", INVALID_VALUE)) != INVALID_VALUE)          arm_swing_gain = value;
    if((value = ini->getd(section, "balance_knee_gain", INVALID_VALUE)) != INVALID_VALUE)       balance_knee_gain = value;
    if((value = ini->getd(section, "balance_ankle_pitch_gain", INVALID_VALUE)) != INVALID_VALUE)balance_ankle_pitch_gain = value;
    if((value = ini->getd(section, "balance_hip_roll_gain", INVALID_VALUE)) != INVALID_VALUE)   balance_hip_roll_gain = value;
    if((value = ini->getd(section, "balance_ankle_roll_gain", INVALID_VALUE)) != INVALID_VALUE) balance_ankle_roll_gain = value;
    int ivalue = INVALID_VALUE;
    if((ivalue = ini->geti(section, "p_gain", INVALID_VALUE)) != INVALID_VALUE)                 p_gain = ivalue;
    if((ivalue = ini->geti(section, "i_gain", INVALID_VALUE)) != INVALID_VALUE)                 i_gain = ivalue;
    if((ivalue = ini->geti(section, "d_gain", INVALID_VALUE)) != INVALID_VALUE)                 d_gain = ivalue;
}
void walking::SaveINISettings(minIni* ini)
{
    SaveINISettings(ini, WALKING_SECTION);
}
void walking::SaveINISettings(minIni* ini, const std::string &section)
{
    ini->put(section,   "x_offset",                 x_offset);
    ini->put(section,   "y_offset",                 y_offset);
    ini->put(section,   "z_offset",                 z_offset);
    ini->put(section,   "roll_offset",              r_offset);
    ini->put(section,   "pitch_offset",             p_offset);
    ini->put(section,   "yaw_offset",               a_offset);
    ini->put(section,   "hip_pitch_offset",         hip_pitch_offset);
    ini->put(section,   "period_time",              period_time);
    ini->put(section,   "dsp_ratio",                dsp_ratio);
    ini->put(section,   "step_forward_back_ratio",  step_fb_ratio);
    ini->put(section,   "foot_height",              z_move_amplitude);
    ini->put(section,   "swing_right_left",         y_swap_amplitude);
    ini->put(section,   "swing_top_down",           z_swap_amplitude);
    ini->put(section,   "pelvis_offset",            pelvis_offset);
    ini->put(section,   "arm_swing_gain",           arm_swing_gain);
    ini->put(section,   "balance_knee_gain",        balance_knee_gain);
    ini->put(section,   "balance_ankle_pitch_gain", balance_ankle_pitch_gain);
    ini->put(section,   "balance_hip_roll_gain",    balance_hip_roll_gain);
    ini->put(section,   "balance_ankle_roll_gain",  balance_ankle_roll_gain);
    ini->put(section,   "p_gain",                   p_gain);
    ini->put(section,   "i_gain",                   i_gain);
    ini->put(section,   "d_gain",                   d_gain);
}*/

void walking::set_angle(int id, double angle) {}


void walking::set_p_gain(int id, int pgain) {}


constexpr double walking::wsin(double time, double period, double period_shift, double mag, double mag_shift) {
    return mag * sin(2 * 3.141592 / period * time - period_shift) + mag_shift;
}


bool walking::compute_ik(double* out, double x, double y, double z, double a, double b, double c) {
    matrix Tad, Tda, Tcd, Tdc, Tac;
    vector vec;
    double _Rac, _Acos, _Atan, _k, _l, _m, _n, _s, _c, _theta;
    double LEG_LENGTH = kinematics::LEG_LENGTH;
    double THIGH_LENGTH = kinematics::THIGH_LENGTH;
    double CALF_LENGTH = kinematics::CALF_LENGTH;
    double ANKLE_LENGTH = kinematics::ANKLE_LENGTH;

    Tad.set_transform(point3d(x, y, z - LEG_LENGTH), vector(a * 180.0 / PI, b * 180.0 / PI, c * 180.0 / PI));

    vec.x = x + Tad.m[2] * ANKLE_LENGTH;
    vec.y = y + Tad.m[6] * ANKLE_LENGTH;
    vec.z = (z - LEG_LENGTH) + Tad.m[10] * ANKLE_LENGTH;

    // Get Knee
    _Rac = vec.length();
    _Acos = acos(
            (_Rac * _Rac - THIGH_LENGTH * THIGH_LENGTH - CALF_LENGTH * CALF_LENGTH) / (2 * THIGH_LENGTH * CALF_LENGTH));
    if (isnan(_Acos) == 1)
        return false;
    *(out + 3) = _Acos;

    // Get Ankle Roll
    Tda = Tad;
    if (!Tda.inverse())
        return false;
    _k = sqrt(Tda.m[7] * Tda.m[7] + Tda.m[11] * Tda.m[11]);
    _l = sqrt(Tda.m[7] * Tda.m[7] + (Tda.m[11] - ANKLE_LENGTH) * (Tda.m[11] - ANKLE_LENGTH));
    _m = (_k * _k - _l * _l - ANKLE_LENGTH * ANKLE_LENGTH) / (2 * _l * ANKLE_LENGTH);
    if (_m > 1.0)
        _m = 1.0;
    else if (_m < -1.0)
        _m = -1.0;
    _Acos = acos(_m);
    if (isnan(_Acos) == 1)
        return false;
    if (Tda.m[7] < 0.0)
        *(out + 5) = -_Acos;
    else
        *(out + 5) = _Acos;

    // Get Hip Yaw
    Tcd.set_transform(point3d(0, 0, -ANKLE_LENGTH), vector(*(out + 5) * 180.0 / PI, 0, 0));
    Tdc = Tcd;
    if (!Tdc.inverse())
        return false;
    Tac = Tad * Tdc;
    _Atan = atan2(-Tac.m[1], Tac.m[5]);
    if (isinf(_Atan)) return false;
    *(out) = _Atan;

    // Get Hip Roll
    _Atan = atan2(Tac.m[9], -Tac.m[1] * sin(*(out)) + Tac.m[5] * cos(*(out)));
    if (isinf(_Atan)) return false;
    *(out + 1) = _Atan;

    // Get Hip Pitch and Ankle Pitch
    _Atan = atan2(Tac.m[2] * cos(*(out)) + Tac.m[6] * sin(*(out)), Tac.m[0] * cos(*(out)) + Tac.m[4] * sin(*(out)));
    if (isinf(_Atan)) return false;
    _theta = _Atan;
    _k = sin(*(out + 3)) * CALF_LENGTH;
    _l = -THIGH_LENGTH - cos(*(out + 3)) * CALF_LENGTH;
    _m = cos(*(out)) * vec.x + sin(*(out)) * vec.y;
    _n = cos(*(out + 1)) * vec.z + sin(*(out)) * sin(*(out + 1)) * vec.x - cos(*(out)) * sin(*(out + 1)) * vec.y;
    _s = (_k * _n + _l * _m) / (_k * _k + _l * _l);
    _c = (_n - _k * _s) / _l;
    _Atan = atan2(_s, _c);
    if (isinf(_Atan)) return false;
    *(out + 2) = _Atan;
    *(out + 4) = _theta - *(out + 3) - *(out + 2);
    return true;
}


void walking::update_param_time() {
    m_period_time = period_time;
    m_dsp_ratio = dsp_ratio;
    m_ssp_ratio = 1 - dsp_ratio;

    m_x_swap_period_time = m_period_time / 2;
    m_x_move_period_time = m_period_time * m_ssp_ratio;
    m_y_swap_period_time = m_period_time;
    m_y_move_period_time = m_period_time * m_ssp_ratio;
    m_z_swap_period_time = m_period_time / 2;
    m_z_move_period_time = m_period_time * m_ssp_ratio / 2;
    m_a_move_period_time = m_period_time * m_ssp_ratio;

    m_ssp_time = m_period_time * m_ssp_ratio;
    m_ssp_time_start_l = (1 - m_ssp_ratio) * m_period_time / 4;
    m_ssp_time_end_l = (1 + m_ssp_ratio) * m_period_time / 4;
    m_ssp_time_start_r = (3 - m_ssp_ratio) * m_period_time / 4;
    m_ssp_time_end_r = (3 + m_ssp_ratio) * m_period_time / 4;

    m_phase_time1 = (m_ssp_time_end_l + m_ssp_time_start_l) / 2;
    m_phase_time2 = (m_ssp_time_start_r + m_ssp_time_end_l) / 2;
    m_phase_time3 = (m_ssp_time_end_r + m_ssp_time_start_r) / 2;

    m_pelvis_offset = pelvis_offset;//*MX28::RATIO_ANGLE2VALUE
    m_pelvis_swing = m_pelvis_offset * 0.35;
    m_arm_swing_gain = arm_swing_gain;
}


void walking::update_param_move() {
    // Forward/Back
    m_x_move_amplitude = x_move_amplitude;
    m_x_swap_amplitude = x_move_amplitude * step_fb_ratio;

    // Right/Left
    m_y_move_amplitude = y_move_amplitude / 2;
    if (m_y_move_amplitude > 0)
        m_y_move_amplitude_shift = m_y_move_amplitude;
    else
        m_y_move_amplitude_shift = -m_y_move_amplitude;
    m_y_swap_amplitude = y_swap_amplitude + m_y_move_amplitude_shift * 0.04;

    m_z_move_amplitude = z_move_amplitude / 2;
    m_z_move_amplitude_shift = m_z_move_amplitude / 2;
    m_z_swap_amplitude = z_swap_amplitude;
    m_z_swap_amplitude_shift = m_z_swap_amplitude;

    // Direction
    if (!a_move_aim_on) {
        m_a_move_amplitude = a_move_amplitude * PI / 180.0 / 2;
        if (m_a_move_amplitude > 0)
            m_a_move_amplitude_shift = m_a_move_amplitude;
        else
            m_a_move_amplitude_shift = -m_a_move_amplitude;
    } else {
        m_a_move_amplitude = -a_move_amplitude * PI / 180.0 / 2;
        if (m_a_move_amplitude > 0)
            m_a_move_amplitude_shift = -m_a_move_amplitude;
        else
            m_a_move_amplitude_shift = m_a_move_amplitude;
    }
}


void walking::update_param_balance() {
    m_x_offset = x_offset;
    m_y_offset = y_offset;
    m_z_offset = z_offset;
    m_r_offset = r_offset * PI / 180.0;
    m_p_offset = p_offset * PI / 180.0;
    m_a_offset = a_offset * PI / 180.0;
    m_hip_pitch_offset = hip_pitch_offset;//*MX28::RATIO_ANGLE2VALUE;//todo
}


void walking::initialize() {
    x_move_amplitude = 0;
    y_move_amplitude = 0;
    a_move_amplitude = 0;

    m_body_swing_y = 0;
    m_body_swing_z = 0;

    m_x_swap_phase_shift = PI;
    m_x_swap_amplitude_shift = 0;
    m_x_move_phase_shift = PI / 2;
    m_x_move_amplitude_shift = 0;
    m_y_swap_phase_shift = 0;
    m_y_swap_amplitude_shift = 0;
    m_y_move_phase_shift = PI / 2;
    m_z_swap_phase_shift = PI * 3 / 2;
    m_z_move_phase_shift = PI / 2;
    m_a_move_phase_shift = PI / 2;

    m_ctrl_running = false;
    m_real_running = false;
    m_time = 0;

    this->update_param_time();
    this->update_param_move();

    this->process();
}


void walking::start() {
    m_ctrl_running = true;
    m_real_running = true;
}


void walking::stop() {
    m_ctrl_running = false;
}


bool walking::is_running() {
    return m_real_running;
}


void walking::process() {
    double x_swap, y_swap, z_swap, a_swap, b_swap, c_swap;
    double x_move_r, y_move_r, z_move_r, a_move_r, b_move_r, c_move_r;
    double x_move_l, y_move_l, z_move_l, a_move_l, b_move_l, c_move_l;
    double pelvis_offset_r, pelvis_offset_l;
    double angle[14], ep[12];
    //                     R_HIP_YAW, R_HIP_ROLL, R_HIP_PITCH, R_KNEE, R_ANKLE_PITCH, R_ANKLE_ROLL, L_HIP_YAW, L_HIP_ROLL, L_HIP_PITCH, L_KNEE, L_ANKLE_PITCH, L_ANKLE_ROLL, R_ARM_SWING, L_ARM_SWING
    double dir[14] = {-1, -1, 1, 1, -1, 1, -1, -1, -1, -1, 1, 1, 1, -1};

    // Update walk parameters
    if (m_time == 0) {
        this->update_param_time();
        if (!m_ctrl_running) {
            if (m_x_move_amplitude == 0 &&
                    m_y_move_amplitude == 0 &&
                    m_a_move_amplitude == 0) {
                m_real_running = false;
            } else {
                x_move_amplitude = 0;
                y_move_amplitude = 0;
                a_move_amplitude = 0;
            }
        }
    } else if (m_time >= (m_phase_time1 - TIME_UNIT / 2) &&
               m_time < (m_phase_time1 + TIME_UNIT / 2)) {
        this->update_param_move();
    } else if (m_time >= (m_phase_time2 - TIME_UNIT / 2) &&
               m_time < (m_phase_time2 + TIME_UNIT / 2)) {
        this->update_param_time();
        m_time = m_phase_time2;
        if (!m_ctrl_running) {
            if (m_x_move_amplitude == 0 &&
                    m_y_move_amplitude == 0 &&
                    m_a_move_amplitude == 0) {
                m_real_running = false;
            } else {
                x_move_amplitude = 0;
                y_move_amplitude = 0;
                a_move_amplitude = 0;
            }
        }
    } else if (m_time >= (m_phase_time3 - TIME_UNIT / 2) &&
               m_time < (m_phase_time3 + TIME_UNIT / 2)) {
        this->update_param_move();
    }
    this->update_param_balance();

    // Compute endpoints
    x_swap = wsin(m_time, m_x_swap_period_time, m_x_swap_phase_shift, m_x_swap_amplitude, m_x_swap_amplitude_shift);
    y_swap = wsin(m_time, m_y_swap_period_time, m_y_swap_phase_shift, m_y_swap_amplitude, m_y_swap_amplitude_shift);
    z_swap = wsin(m_time, m_z_swap_period_time, m_z_swap_phase_shift, m_z_swap_amplitude, m_z_swap_amplitude_shift);
    a_swap = 0;
    b_swap = 0;
    c_swap = 0;

    if (m_time <= m_ssp_time_start_l) {
        x_move_l = wsin(m_ssp_time_start_l,
                        m_x_move_period_time,
                        m_x_move_phase_shift + 2 * PI / m_x_move_period_time * m_ssp_time_start_l,
                        m_x_move_amplitude,
                        m_x_move_amplitude_shift);
        y_move_l = wsin(m_ssp_time_start_l,
                        m_y_move_period_time,
                        m_y_move_phase_shift + 2 * PI / m_y_move_period_time * m_ssp_time_start_l,
                        m_y_move_amplitude,
                        m_y_move_amplitude_shift);
        z_move_l = wsin(m_ssp_time_start_l,
                        m_z_move_period_time,
                        m_z_move_phase_shift + 2 * PI / m_z_move_period_time * m_ssp_time_start_l,
                        m_z_move_amplitude,
                        m_z_move_amplitude_shift);
        c_move_l = wsin(m_ssp_time_start_l,
                        m_a_move_period_time,
                        m_a_move_phase_shift + 2 * PI / m_a_move_period_time * m_ssp_time_start_l,
                        m_a_move_amplitude,
                        m_a_move_amplitude_shift);
        x_move_r = wsin(m_ssp_time_start_l,
                        m_x_move_period_time,
                        m_x_move_phase_shift + 2 * PI / m_x_move_period_time * m_ssp_time_start_l,
                        -m_x_move_amplitude,
                        -m_x_move_amplitude_shift);
        y_move_r = wsin(m_ssp_time_start_l,
                        m_y_move_period_time,
                        m_y_move_phase_shift + 2 * PI / m_y_move_period_time * m_ssp_time_start_l,
                        -m_y_move_amplitude,
                        -m_y_move_amplitude_shift);
        z_move_r = wsin(m_ssp_time_start_r,
                        m_z_move_period_time,
                        m_z_move_phase_shift + 2 * PI / m_z_move_period_time * m_ssp_time_start_r,
                        m_z_move_amplitude,
                        m_z_move_amplitude_shift);
        c_move_r = wsin(m_ssp_time_start_l, m_a_move_period_time,
                        m_a_move_phase_shift + 2 * PI / m_a_move_period_time * m_ssp_time_start_l, -m_a_move_amplitude,
                        -m_a_move_amplitude_shift);
        pelvis_offset_l = 0;
        pelvis_offset_r = 0;
    } else if (m_time <= m_ssp_time_end_l) {
        x_move_l = wsin(m_time,
                        m_x_move_period_time,
                        m_x_move_phase_shift + 2 * PI / m_x_move_period_time * m_ssp_time_start_l,
                        m_x_move_amplitude,
                        m_x_move_amplitude_shift);
        y_move_l = wsin(m_time,
                        m_y_move_period_time,
                        m_y_move_phase_shift + 2 * PI / m_y_move_period_time * m_ssp_time_start_l,
                        m_y_move_amplitude,
                        m_y_move_amplitude_shift);
        z_move_l = wsin(m_time,
                        m_z_move_period_time,
                        m_z_move_phase_shift + 2 * PI / m_z_move_period_time * m_ssp_time_start_l,
                        m_z_move_amplitude,
                        m_z_move_amplitude_shift);
        c_move_l = wsin(m_time,
                        m_a_move_period_time,
                        m_a_move_phase_shift + 2 * PI / m_a_move_period_time * m_ssp_time_start_l,
                        m_a_move_amplitude,
                        m_a_move_amplitude_shift);
        x_move_r = wsin(m_time,
                        m_x_move_period_time,
                        m_x_move_phase_shift + 2 * PI / m_x_move_period_time * m_ssp_time_start_l,
                        -m_x_move_amplitude,
                        -m_x_move_amplitude_shift);
        y_move_r = wsin(m_time,
                        m_y_move_period_time,
                        m_y_move_phase_shift + 2 * PI / m_y_move_period_time * m_ssp_time_start_l,
                        -m_y_move_amplitude,
                        -m_y_move_amplitude_shift);
        z_move_r = wsin(m_ssp_time_start_r,
                        m_z_move_period_time,
                        m_z_move_phase_shift + 2 * PI / m_z_move_period_time * m_ssp_time_start_r,
                        m_z_move_amplitude,
                        m_z_move_amplitude_shift);
        c_move_r = wsin(m_time,
                        m_a_move_period_time,
                        m_a_move_phase_shift + 2 * PI / m_a_move_period_time * m_ssp_time_start_l,
                        -m_a_move_amplitude,
                        -m_a_move_amplitude_shift);
        pelvis_offset_l = wsin(m_time,
                               m_z_move_period_time,
                               m_z_move_phase_shift + 2 * PI / m_z_move_period_time * m_ssp_time_start_l,
                               m_pelvis_swing / 2,
                               m_pelvis_swing / 2);
        pelvis_offset_r = wsin(m_time,
                               m_z_move_period_time,
                               m_z_move_phase_shift + 2 * PI / m_z_move_period_time * m_ssp_time_start_l,
                               -m_pelvis_offset / 2,
                               -m_pelvis_offset / 2);
    } else if (m_time <= m_ssp_time_start_r) {
        x_move_l = wsin(m_ssp_time_end_l,
                        m_x_move_period_time,
                        m_x_move_phase_shift + 2 * PI / m_x_move_period_time * m_ssp_time_start_l,
                        m_x_move_amplitude,
                        m_x_move_amplitude_shift);
        y_move_l = wsin(m_ssp_time_end_l,
                        m_y_move_period_time,
                        m_y_move_phase_shift + 2 * PI / m_y_move_period_time * m_ssp_time_start_l,
                        m_y_move_amplitude,
                        m_y_move_amplitude_shift);
        z_move_l = wsin(m_ssp_time_end_l,
                        m_z_move_period_time,
                        m_z_move_phase_shift + 2 * PI / m_z_move_period_time * m_ssp_time_start_l,
                        m_z_move_amplitude,
                        m_z_move_amplitude_shift);
        c_move_l = wsin(m_ssp_time_end_l,
                        m_a_move_period_time,
                        m_a_move_phase_shift + 2 * PI / m_a_move_period_time * m_ssp_time_start_l,
                        m_a_move_amplitude,
                        m_a_move_amplitude_shift);
        x_move_r = wsin(m_ssp_time_end_l,
                        m_x_move_period_time,
                        m_x_move_phase_shift + 2 * PI / m_x_move_period_time * m_ssp_time_start_l,
                        -m_x_move_amplitude,
                        -m_x_move_amplitude_shift);
        y_move_r = wsin(m_ssp_time_end_l,
                        m_y_move_period_time,
                        m_y_move_phase_shift + 2 * PI / m_y_move_period_time * m_ssp_time_start_l,
                        -m_y_move_amplitude,
                        -m_y_move_amplitude_shift);
        z_move_r = wsin(m_ssp_time_start_r,
                        m_z_move_period_time,
                        m_z_move_phase_shift + 2 * PI / m_z_move_period_time * m_ssp_time_start_r,
                        m_z_move_amplitude,
                        m_z_move_amplitude_shift);
        c_move_r = wsin(m_ssp_time_end_l,
                        m_a_move_period_time,
                        m_a_move_phase_shift + 2 * PI / m_a_move_period_time * m_ssp_time_start_l,
                        -m_a_move_amplitude,
                        -m_a_move_amplitude_shift);
        pelvis_offset_l = 0;
        pelvis_offset_r = 0;
    } else if (m_time <= m_ssp_time_end_r) {
        x_move_l = wsin(m_time,
                        m_x_move_period_time,
                        m_x_move_phase_shift + 2 * PI / m_x_move_period_time * m_ssp_time_start_r + PI,
                        m_x_move_amplitude,
                        m_x_move_amplitude_shift);
        y_move_l = wsin(m_time,
                        m_y_move_period_time,
                        m_y_move_phase_shift + 2 * PI / m_y_move_period_time * m_ssp_time_start_r + PI,
                        m_y_move_amplitude,
                        m_y_move_amplitude_shift);
        z_move_l = wsin(m_ssp_time_end_l,
                        m_z_move_period_time,
                        m_z_move_phase_shift + 2 * PI / m_z_move_period_time * m_ssp_time_start_l,
                        m_z_move_amplitude,
                        m_z_move_amplitude_shift);
        c_move_l = wsin(m_time,
                        m_a_move_period_time,
                        m_a_move_phase_shift + 2 * PI / m_a_move_period_time * m_ssp_time_start_r + PI,
                        m_a_move_amplitude,
                        m_a_move_amplitude_shift);
        x_move_r = wsin(m_time,
                        m_x_move_period_time,
                        m_x_move_phase_shift + 2 * PI / m_x_move_period_time * m_ssp_time_start_r + PI,
                        -m_x_move_amplitude,
                        -m_x_move_amplitude_shift);
        y_move_r = wsin(m_time,
                        m_y_move_period_time,
                        m_y_move_phase_shift + 2 * PI / m_y_move_period_time * m_ssp_time_start_r + PI,
                        -m_y_move_amplitude,
                        -m_y_move_amplitude_shift);
        z_move_r = wsin(m_time,
                        m_z_move_period_time,
                        m_z_move_phase_shift + 2 * PI / m_z_move_period_time * m_ssp_time_start_r,
                        m_z_move_amplitude,
                        m_z_move_amplitude_shift);
        c_move_r = wsin(m_time,
                        m_a_move_period_time,
                        m_a_move_phase_shift + 2 * PI / m_a_move_period_time * m_ssp_time_start_r + PI,
                        -m_a_move_amplitude,
                        -m_a_move_amplitude_shift);
        pelvis_offset_l = wsin(m_time,
                               m_z_move_period_time,
                               m_z_move_phase_shift + 2 * PI / m_z_move_period_time * m_ssp_time_start_r,
                               m_pelvis_offset / 2,
                               m_pelvis_offset / 2);
        pelvis_offset_r = wsin(m_time,
                               m_z_move_period_time,
                               m_z_move_phase_shift + 2 * PI / m_z_move_period_time * m_ssp_time_start_r,
                               -m_pelvis_swing / 2,
                               -m_pelvis_swing / 2);
    } else {
        x_move_l = wsin(m_ssp_time_end_r,
                        m_x_move_period_time,
                        m_x_move_phase_shift + 2 * PI / m_x_move_period_time * m_ssp_time_start_r + PI,
                        m_x_move_amplitude,
                        m_x_move_amplitude_shift);
        y_move_l = wsin(m_ssp_time_end_r,
                        m_y_move_period_time,
                        m_y_move_phase_shift + 2 * PI / m_y_move_period_time * m_ssp_time_start_r + PI,
                        m_y_move_amplitude,
                        m_y_move_amplitude_shift);
        z_move_l = wsin(m_ssp_time_end_l,
                        m_z_move_period_time,
                        m_z_move_phase_shift + 2 * PI / m_z_move_period_time * m_ssp_time_start_l,
                        m_z_move_amplitude,
                        m_z_move_amplitude_shift);
        c_move_l = wsin(m_ssp_time_end_r,
                        m_a_move_period_time,
                        m_a_move_phase_shift + 2 * PI / m_a_move_period_time * m_ssp_time_start_r + PI,
                        m_a_move_amplitude,
                        m_a_move_amplitude_shift);
        x_move_r = wsin(m_ssp_time_end_r,
                        m_x_move_period_time,
                        m_x_move_phase_shift + 2 * PI / m_x_move_period_time * m_ssp_time_start_r + PI,
                        -m_x_move_amplitude,
                        -m_x_move_amplitude_shift);
        y_move_r = wsin(m_ssp_time_end_r,
                        m_y_move_period_time,
                        m_y_move_phase_shift + 2 * PI / m_y_move_period_time * m_ssp_time_start_r + PI,
                        -m_y_move_amplitude,
                        -m_y_move_amplitude_shift);
        z_move_r = wsin(m_ssp_time_end_r,
                        m_z_move_period_time,
                        m_z_move_phase_shift + 2 * PI / m_z_move_period_time * m_ssp_time_start_r,
                        m_z_move_amplitude,
                        m_z_move_amplitude_shift);
        c_move_r = wsin(m_ssp_time_end_r,
                        m_a_move_period_time,
                        m_a_move_phase_shift + 2 * PI / m_a_move_period_time * m_ssp_time_start_r + PI,
                        -m_a_move_amplitude,
                        -m_a_move_amplitude_shift);
        pelvis_offset_l = 0;
        pelvis_offset_r = 0;
    }

    a_move_l = 0;
    b_move_l = 0;
    a_move_r = 0;
    b_move_r = 0;

    ep[R_X] = x_swap + x_move_r + m_x_offset;
    ep[R_Y] = y_swap + y_move_r - m_y_offset / 2;
    ep[R_Z] = z_swap + z_move_r + m_z_offset;
    ep[R_ROLL] = a_swap + a_move_r - m_r_offset / 2;
    ep[R_PITCH] = b_swap + b_move_r + m_p_offset;
    ep[R_YAW] = c_swap + c_move_r - m_a_offset / 2;
    ep[L_X] = x_swap + x_move_l + m_x_offset;
    ep[L_Y] = y_swap + y_move_l + m_y_offset / 2;
    ep[L_Z] = z_swap + z_move_l + m_z_offset;
    ep[L_ROLL] = a_swap + a_move_l + m_r_offset / 2;
    ep[L_PITCH] = b_swap + b_move_l + m_p_offset;
    ep[L_YAW] = c_swap + c_move_l + m_a_offset / 2;

    // Compute body swing
    if (m_time <= m_ssp_time_end_l) {
        m_body_swing_y = -ep[L_HIP_ROLL]; //[L_Y]
        m_body_swing_z = ep[L_HIP_PITCH]; //[L_Z]
    } else {
        m_body_swing_y = -ep[R_HIP_ROLL];//[R_Y]
        m_body_swing_z = ep[R_HIP_PITCH];//[R_Z]
    }
    m_body_swing_z -= kinematics::LEG_LENGTH;

    // Compute arm swing
    if (m_x_move_amplitude == 0) {
        angle[R_ARM_SWING] = 0;
        angle[L_ARM_SWING] = 0;
    } else {
        angle[R_ARM_SWING] = wsin(m_time, m_period_time, PI * 1.5, -m_x_move_amplitude * m_arm_swing_gain, 0);
        angle[L_ARM_SWING] = wsin(m_time, m_period_time, PI * 1.5, m_x_move_amplitude * m_arm_swing_gain, 0);
    }

    if (m_real_running) {
        m_time += TIME_UNIT;
        if (m_time >= m_period_time)
            m_time = 0;
    }	

    // Compute angles
    if (!this->compute_ik(&angle[R_HIP_YAW], ep[R_X], ep[R_Y], ep[R_Z], ep[R_ROLL], ep[R_PITCH], ep[R_YAW]) ||
            !this->compute_ik(&angle[L_HIP_YAW], ep[L_X], ep[L_Y], ep[L_Z], ep[L_ROLL], ep[L_PITCH], ep[L_YAW])) {
        return; // Do not use angle;
    }

    for (int i = 0; i < 12; i++)
        angle[i] *= 180.0 / PI;

    // Compute motor value	
/*pelvis_offset_r pelvis_offset_l показывают градус поворота hip_pitch_offset показывает градус наклона торса*/
double sin_r = sin(ep[R_YAW]);
double sin_l = sin(ep[L_YAW]);
double cos_r = cos(ep[R_YAW]);
double cos_l = cos(ep[L_YAW]);

    angle[R_HIP_ROLL] += hip_pitch_offset * sin_r + pelvis_offset_r * cos_r;
    angle[L_HIP_ROLL] += hip_pitch_offset * sin_l + pelvis_offset_l * cos_l;
    angle[R_HIP_PITCH] -= hip_pitch_offset * cos_r + pelvis_offset_r * sin_r;
    angle[L_HIP_PITCH] -= hip_pitch_offset * cos_l + pelvis_offset_l * sin_l;
/*   
//pitch(вокруг оси у) Yaw (z)
    angle[R_HIP_ROLL] += m_body_swing_z*sin(ep[R_YAW])+m_body_swing_y*cos(ep[R_YAW]);
    angle[L_HIP_ROLL] += m_body_swing_z*sin(ep[L_YAW])+m_body_swing_y*cos(ep[L_YAW]);
    angle[R_HIP_PITCH] -= m_body_swing_z*cos(ep[R_YAW])+m_body_swing_y*sin(ep[R_YAW]);
    angle[L_HIP_PITCH] -= m_body_swing_z*cos(ep[L_YAW])+m_body_swing_y*sin(ep[L_YAW]);

*/


    
    angle[R_ARM_SWING] -= 48.345;
    angle[L_ARM_SWING] += 41.313;

    for (int i = 0; i < NUMBER_OF_JOINTS; i++) {
        angle[i] = angle[i] * dir[i];
    }

    // adjust balance offset
    if (m_balance_enable) {
//         TODO Read from message
        double rl_gyro_err = 0;// RL_GYRO;
        double fb_gyro_err = 0;// FB_GYRO;

        // TODO TUNE IT!
/*
        angle[R_HIP_ROLL] += dir[R_HIP_ROLL] * rl_gyro_err * balance_hip_roll_gain * 4;//1
        angle[L_HIP_ROLL] += dir[L_HIP_ROLL] * rl_gyro_err * balance_hip_roll_gain * 4;//7
*/
        angle[R_HIP_ROLL] += dir[R_HIP_ROLL] * (rl_gyro_err * balance_hip_roll_gain * cos_r + fb_gyro_err * FB_BALANCE_HIP_PITCH_GAIN * sin_r )* 4;//1
        angle[L_HIP_ROLL] += dir[L_HIP_ROLL] * (rl_gyro_err * balance_hip_roll_gain * cos_l + fb_gyro_err * FB_BALANCE_HIP_PITCH_GAIN * sin_r ) * 4;//7



        
//added
        angle[R_HIP_PITCH] += dir[R_HIP_PITCH] * (rl_gyro_err * balance_hip_roll_gain * sin_r + fb_gyro_err * FB_BALANCE_HIP_PITCH_GAIN * cos_r ) * 4;//2
        angle[L_HIP_PITCH] += dir[L_HIP_PITCH] * (rl_gyro_err * balance_hip_roll_gain * sin_l + fb_gyro_err * FB_BALANCE_HIP_PITCH_GAIN * cos_l ) * 4;//8

        angle[R_KNEE] -= dir[R_KNEE] * fb_gyro_err * balance_knee_gain * 4;//3
        angle[L_KNEE] -= dir[L_KNEE] * fb_gyro_err * balance_knee_gain * 4;//9
        angle[R_ANKLE_PITCH] -= dir[R_ANKLE_PITCH] * fb_gyro_err * balance_ankle_pitch_gain * 4;//4
        angle[L_ANKLE_PITCH] -= dir[L_ANKLE_PITCH] * fb_gyro_err * balance_ankle_pitch_gain * 4;//10
        angle[R_ANKLE_ROLL] -= dir[R_ANKLE_ROLL] * rl_gyro_err * balance_ankle_roll_gain * 4;//5
        angle[L_ANKLE_ROLL] -= dir[L_ANKLE_ROLL] * rl_gyro_err * balance_ankle_roll_gain * 4;//11
    }

    this->set_angle(joint::ID_R_HIP_YAW, angle[0]);
    this->set_angle(joint::ID_R_HIP_ROLL, angle[1]);
    this->set_angle(joint::ID_R_HIP_PITCH, angle[2]);
    this->set_angle(joint::ID_R_KNEE, angle[3]);
    this->set_angle(joint::ID_R_ANKLE_PITCH, angle[4]);
    this->set_angle(joint::ID_R_ANKLE_ROLL, angle[5]);
    this->set_angle(joint::ID_L_HIP_YAW, angle[6]);
    this->set_angle(joint::ID_L_HIP_ROLL, angle[7]);
    this->set_angle(joint::ID_L_HIP_PITCH, angle[8]);
    this->set_angle(joint::ID_L_KNEE, angle[9]);
    this->set_angle(joint::ID_L_ANKLE_PITCH, angle[10]);
    this->set_angle(joint::ID_L_ANKLE_ROLL, angle[11]);
    this->set_angle(joint::ID_R_SHOULDER_PITCH, angle[12]);
    this->set_angle(joint::ID_L_SHOULDER_PITCH, angle[13]);

    //this->set_angle(joint::ID_HEAD_PAN, a_move_amplitude); // O_O

    /* todo
     for(int id = JointData::ID_R_HIP_YAW; id <= JointData::ID_L_ANKLE_ROLL; id++)
    {
        m_Joint.SetPGain(id, p_gain);
        m_Joint.SetIGain(id, i_gain);
        m_Joint.SetDGain(id, d_gain);
    }
    */
}

    Contact GitHub API Training Shop Blog About 

    © 2017 GitHub, Inc. Terms Privacy Security Status Help 



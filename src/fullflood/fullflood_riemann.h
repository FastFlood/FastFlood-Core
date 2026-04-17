#ifndef FULLFLOOD_RIEMANN_H
#define FULLFLOOD_RIEMANN_H


static LSMVector3 F_HLL2FL(float h_L,float u_L,float v_L,float h_R,float u_R,float v_R)
{
        LSMVector3 ret;

        float grav_h_L = GRAV*h_L;
        float grav_h_R = GRAV*h_R;
        float sqrt_grav_h_L = std::sqrt((float)(grav_h_L));  // wave velocity
        float sqrt_grav_h_R = std::sqrt((float)(grav_h_R));
        float q_R = u_R*h_R;
        float q_L = u_L*h_L;

        float sc1 = std::min((float)(u_L - sqrt_grav_h_L),(float)(u_R - sqrt_grav_h_R));
        float sc2 = std::max((float)(u_L + sqrt_grav_h_L),(float)(u_R + sqrt_grav_h_R));
        float tmp = 1.0/std::max((float)0.1,(float)(sc2-sc1));
        float t1 = (std::min((float)(sc2),(float)(0.)) - std::min((float)(sc1),(float)(0.)))*tmp;
        float t2 = 1.0 - t1;
        float t3_1 = std::fabs((float)(sc1));
        float t3_2 = std::fabs((float)(sc2));
        float t3 =(sc2* t3_1 - sc1* t3_2)*0.5*tmp;

        if(h_L > 0.0 || h_R > 0.0 )
        {

            //FF::cout << h_L << " " << h_R << " "<< t1 << " "  << t2 << std::endl;
        }
        ret.x = t1*q_R + t2*q_L - t3*(h_R - h_L);
        ret.y = t1*(q_R*u_R + grav_h_R*h_R*0.5) + t2*(q_L*u_L + grav_h_L*h_L*0.5) - t3*(q_R - q_L);
        ret.z = t1*q_R*v_R + t2*q_L*v_L - t3*(h_R*v_R - h_L*v_L);

        return ret;
}




#endif // FULLFLOOD_RIEMANN_H

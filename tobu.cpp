#include "pwm_uart.hpp"
#include "ekf_sensor.hpp"


//e エレベータ a エルロン r ラダー t スロットル
float h=0.01;
//float kpe=0.053051647,kpa=0.053051647,kpr=0.159154943;
float kpe=0.50,kpa=0.50,kpr=0.75;
float shigma_e,shiguma_a,shiguma_r;
float ref_e,ref_r,ref_a,err_e,err_a,err_r;
float sk_e,sk_a,sk_r;
float olderr_e,olderr_a,olderr_r;
float dk_a,dk_e,dk_r;
float se,sa,sr;
float Ti=10000;
float Td=0;
/* Private macro -------------------------------------------------------------*/

int main(void)
{
  Matrix<float, 7 ,1> xp = MatrixXf::Zero(7,1);
  Matrix<float, 7 ,1> xe = MatrixXf::Zero(7,1);
  Matrix<float, 7 ,1> x_sim = MatrixXf::Zero(7,1);
  Matrix<float, 7 ,7> P = MatrixXf::Identity(7,7);
  Matrix<float, 6 ,1> z = MatrixXf::Zero(6,1);
  Matrix<float, 6 ,1> z_sim = MatrixXf::Zero(6,1);
  Matrix<float, 6 ,1> z_noise = MatrixXf::Zero(6,1);
  Matrix<float, 3, 1> omega_m = MatrixXf::Zero(3, 1);
  Matrix<float, 3, 1> omega_sim;
  Matrix<float, 3, 1> domega;
  Matrix<float, 3, 1> domega_sim;
  Matrix<float, 3, 3> Q = MatrixXf::Identity(3, 3)*1;
  Matrix<float, 6, 6> R = MatrixXf::Identity(6, 6)*1;
  Matrix<float, 7 ,3> G;
  Matrix<float, 3 ,1> beta;
  float t=0.0,dt=0.01;
  uint64_t s_time=0,e_time=0,d_time=0;
  double phl,theta,psi;
  short i,waittime=5;
  float p_com[10]={0.1*PI, 0, 0.01*PI, 0, -0.01*PI, 0, 0.02*PI, 0, -0.05*PI, 0};
  float q_com[10]={0.1*PI, 0, 0      , 0, -0.01*PI, 0, 0.02*PI, 0, -0.02*PI, 0};
  float r_com[10]={0.1*PI, 0,-0.02*PI, 0,      -PI, 0, 0.02*PI, 0,  0.02*PI, 0};
  float endtime=10000.0;
  float control_period=5.0;
  float control_time=5.0;
  int counter=0;
  int sample=1;
  int control_counter=0;
  int control_counter_max=0;
  double pi=3.14159265358;
  float ax,ay,az,wp,wq,wr,mx,my,mz,wqa=0,wpa=0,wra=0,Wqa,Wpa,Wra; 
  float p11=-0.60526553,p12=0.79021444,p13=-0.09599364,p21=0.78892428,p22=0.61155945,p23=0.05994594,p31=-0.10607597,p32=0.0394485,p33=0.99357521;
  float r1=-4.96008102e-6,r2=-4.60371715e-6,r3=-4.17649591e-6;
  float w,f=1;
  float dmx,dmy,dmz,dxx,dxy,dxz;
  float ddxx, ddxy,ddxz;
  float dddxx,dddxy,dddxz;
  float ddddxx,ddddxy,ddddxz;
 //pwm_uart
  float duty_rr,duty_rl,duty_fl,duty_fr; 
  float Duty_rr,Duty_rl,Duty_fl,Duty_fr;
  float seigyo=0.5;  
  const uint LED_PIN = 25;          //LED_PIN=0
  
  gpio_init(LED_PIN);                       //gpioを使えるようにする
  gpio_set_dir(LED_PIN, GPIO_OUT);     
  xe << 1.0, 0.0, 0.0, 0.0, -0.078, 0.0016, 0.00063;
  xp =xe;

  G <<  0.0,0.0,0.0, 
        0.0,0.0,0.0, 
        0.0,0.0,0.0, 
        0.0,0.0,0.0, 
        1.0,0.0,0.0, 
        0.0,1.0,0.0, 
        0.0,0.0,1.0;

  beta << 0.003, 0.003, 0.003;

  P <<  1,0,0,0,0,0,0,  
        0,1,0,0,0,0,0,
        0,0,1,0,0,0,0,  
        0,0,0,1,0,0,0, 
        0,0,0,0,1,0,0,  
        0,0,0,0,0,1,0,  
        0,0,0,0,0,0,1;

  Q << 7.34944e-6,0,0,
       0,6.861e-6,0,
       0,0,5.195e-6;

  R << 3.608e-6,0,0,0,0,0,
       0,6.261e-6,0,0,0,0,
       0,0,1.889e-5,0,0,0,
       0,0,0,1.0,0,0,
       0,0,0,0,1.0,0,
       0,0,0,0,0,1.0;


  //gpio_put(LED_PIN, 1);
  stdio_init_all();
  imu_mag_init();
  pwm_settei();
  serial_settei();
  
 /* for (i=0;i<waittime;i++)
  {
    printf("#Please wait %d[s] ! \n",waittime-i);
    sleep_ms(1000);
  }
  printf("#Start Kalman Filter\n");
*/

  while(f<=400){
    imu_mag_data_read();
    wp=    angular_rate_mdps[0]*0.001*0.017453292;
    wq=    angular_rate_mdps[1]*0.001*0.017453292;
    wr=   -angular_rate_mdps[2]*0.001*0.017453292;
    wqa=wq+wqa;
    wpa=wp+wpa;
    wra=wr+wra;
    f=f+1;
    sleep_us(1250);
  }
  Wqa=wqa/400;
  Wpa=wpa/400;
  Wra=wra/400;

  //printf("finish average\n");
  gpio_put(LED_PIN, 1);
  while(1)/*t<endtime*/
  {
   // s_time=time_us_64();
    //sleep_ms(10);
    //Control
 /*   if(t>control_time)
    {
      control_time = control_time + control_period;
      control_counter++;
      if(control_counter>control_counter_max)control_counter=0;
    }*/
    //gpio_put(LED_PIN, 1);
    imu_mag_data_read();
    ax=   -acceleration_mg[0]*0.001*GRAV;
    ay=   -acceleration_mg[1]*0.001*GRAV;
    az=    acceleration_mg[2]*0.001*GRAV;
    wp=    angular_rate_mdps[0]*0.001*0.017453292-Wpa;
    wq=    angular_rate_mdps[1]*0.001*0.017453292-Wqa;
    wr=   -angular_rate_mdps[2]*0.001*0.017453292-Wra;
    dmx=  -(magnetic_field_mgauss[0]-310);
    dmy=   (magnetic_field_mgauss[1]-10);
    dmz=  -(magnetic_field_mgauss[2]);
   
   
     
#if 0
    // elevator();    
    //校正作業
    dxx=p11*dmx+p21*dmy+p31*dmz;
    dxy=p12*dmx+p22*dmy+p
      32*dmz;
    dxz=p13*dmx+p23*dmy+p33*dmz;

    ddxx=dxx-22.760831749415342;
    ddxy=dxy-19.734355196006327;
    ddxz=dxz-141.33565570453044;


    w=-1.087745370038146;
  
    dddxx=ddxx*0.0020572671658147883;
    dddxy=ddxy*0.0021354074993493823;
    dddxz=ddxz*0.0019594870993397107;

    mx=p11*dddxx+p22*dddxy+p13*dddxz;
    my=p21*dddxx+p22*dddxy+p23*dddxz;
    mz=p31*dddxx+p32*dddxy+p33*dddxz;
 
    omega_m <<wp,wq,wr;
    z       <<ax,ay,az,mx,my,mz;//ここに入れる
    //--Begin Extended Kalman Filter--
    ekf(xp, xe, P, z, omega_m, Q, R, G*dt, beta, dt);
    if(counter%sample==0)
    {
	              imu_mag_data_read();
                phl=Phl(xe);
                theta=Theta(xe);
                psi=Psi(xe);
        
    
   }
#endif
#if 1    
    //姿勢安定化
    //最大角速度 e,a 6π,r 2π
    //最大角度30°
    ref_e=Data2*18.84955592*0.5;
    ref_a=Data4*18.84955592*0.5;
    ref_r=Data1*6.283185307*0.5;

    //エレベータq
    olderr_e=err_e;
    err_e=(ref_e - wq);
    if (sk_e<=30000){
      sk_e=sk_e+olderr_e;
    }
    else if(-30000<=sk_e){
      sk_e=sk_e+olderr_e;
    }
    dk_e=(err_e-olderr_e)*400;
    se=kpe*(err_e+1/Ti*sk_e+Td*dk_e);

    //エルロンp
 
    olderr_a=err_a;
    err_a=(ref_a - wp);
    if (sk_a<=30000){
      sk_a=sk_a+olderr_a;
    }
    else if(-30000<=sk_a){
      sk_a=sk_a+olderr_a;
    }
    dk_a=(err_a-olderr_a)*400;
    sa=kpa*(err_a+1/Ti*sk_a+Td*dk_a);

    //ラダーr
        
    olderr_r=err_r;
    err_r=(ref_r - wr);
    if (sk_r<=30000){
      sk_r=sk_r+olderr_r;
    }
    else if(-30000<=sk_r){
      sk_r=sk_r+olderr_r;
    }
    dk_r=(err_r-olderr_r)*400;
    sr=kpr*(err_r+1/Ti*sk_r+Td*dk_r);


    Duty_fr=Data3+(se-sa+sr)*0.25;
    Duty_fl=Data3+(se+sa-sr)*0.25;
    Duty_rr=Data3+(-se-sa-sr)*0.25;         
    Duty_rl=Data3+(-se+sa+sr)*0.25;
#endif
#if 0   
    Duty_fr=Data3;
    Duty_fl=Data3;
    Duty_rr=Data3;         
    Duty_rl=Data3;
#endif
    tight_loop_contents();




    duty_rr=(float)(DUTYMAX-DUTYMIN)*Duty_rr+DUTYMIN+53;
    duty_fr=(float)(DUTYMAX-DUTYMIN)*Duty_fr+DUTYMIN;
    duty_rl=(float)(DUTYMAX-DUTYMIN)*Duty_rl+DUTYMIN;
    duty_fl=(float)(DUTYMAX-DUTYMIN)*Duty_fl+DUTYMIN+53;
    //if (duty_rr<DUTYMIN+100)duty_rr=DUTYMIN+5.0;
    //if (duty_rl<DUTYMIN+100)duty_rl=DUTYMIN+5.0;
    //if (duty_fr<DUTYMIN+100)duty_fr=DUTYMIN+5.0;
    //if (duty_fl<DUTYMIN+100)duty_fl=DUTYMIN+5.0;
    
    
    
    
    if (duty_rr>DUTYMAX-50.0)duty_rr=DUTYMAX-50.0;
    if (duty_rr<DUTYMIN+5.0)duty_rr=DUTYMIN+5.0;
    if (duty_fr>DUTYMAX-50.0)duty_fr=DUTYMAX-50.0;
    if (duty_fr<DUTYMIN+5.0)duty_fr=DUTYMIN+5.0;
    if (duty_rl>DUTYMAX-50.0)duty_rl=DUTYMAX-50.0;
    if (duty_rl<DUTYMIN+5.0)duty_rl=DUTYMIN+5.0;
    if (duty_fl>DUTYMAX-50.0)duty_fl=DUTYMAX-50.0;
    if (duty_fl<DUTYMIN+5.0)duty_fl=DUTYMIN+5.0;
   
    if(Data3<0.05)
    {
      pwm_set_chan_level(slice_num[0], PWM_CHAN_A, DUTYMIN);
      pwm_set_chan_level(slice_num[0], PWM_CHAN_B, DUTYMIN);
      pwm_set_chan_level(slice_num[1], PWM_CHAN_A, DUTYMIN);
      pwm_set_chan_level(slice_num[1], PWM_CHAN_B, DUTYMIN);
      sk_r=0;
      sk_a=0;
      sk_e=0;


    }
    else
    {
      pwm_set_chan_level(slice_num[0], PWM_CHAN_A, duty_rr);
      pwm_set_chan_level(slice_num[0], PWM_CHAN_B, duty_fr);
      pwm_set_chan_level(slice_num[1], PWM_CHAN_A, duty_rl);
      pwm_set_chan_level(slice_num[1], PWM_CHAN_B, duty_fl);
    }
    //printf("%04f %04f %04f %04f %04f %04f \n",Olddata[0],Olddata[1],Olddata[2],Olddata[3],Olddata[4],Olddata[5]);    
    sleep_us(1750);
    //e_time=time_us_64();
    //d_time=e_time-s_time;
    //printf("右前%04f   左前%04f   右後ろ%04f   左後ろ%04f\n",duty_fr,duty_fl,duty_rr,duty_rl);
    //printf("%04f %04f %04f %04f \n",t,wp,wq,wr);
    //t=t+1;
  }
}


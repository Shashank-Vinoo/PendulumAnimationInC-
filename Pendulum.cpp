#include <iostream>
#include <cmath>

//creating the state of the double pendulum
struct State{
    double theta1;
    double theta2;
    double omega1;
    double omega2;
};
//creating the params of the double pendulum
struct Params{
    double l1;
    double l2;
    double m1;
    double m2;
    double b1;
    double b2;
    double g;
};
// the four equations of motion — return the rate of change of each state variable
double dtheta1(State s, Params p)  {
    return s.omega1;
}
double dtheta2(State s, Params p){
    return s.omega2;
}
double domega1(State s, Params p) {
    double dtheta = s.theta1 - s.theta2;
    double num = -p.m2 * p.g * std::sin(s.theta1 - 2.0*s.theta2)
                 - 2.0 * std::sin(dtheta) * p.m2 * (p.l2 * s.omega2*s.omega2 + p.l1 * s.omega1*s.omega1 * std::cos(dtheta))
                 - (2.0*p.m1 + p.m2) * p.g * std::sin(s.theta1);
    double den = p.l1 * (2.0*p.m1 + p.m2 - p.m2 * std::cos(2.0*dtheta));
    return num / den - p.b1 * s.omega1;
}
double domega2(State s, Params p) {
    double dtheta = s.theta1 - s.theta2;
    double num = 2.0 * std::sin(dtheta) * ((p.m1 + p.m2) * (p.l1 * s.omega1*s.omega1 + p.g * std::cos(s.theta1))
                 + p.m2 * p.l2 * s.omega2*s.omega2 * std::cos(dtheta));
    double den = p.l2 * (2.0*p.m1 + p.m2 - p.m2 * std::cos(2.0*dtheta));
    return num / den - p.b2 * s.omega2;
}
// advances the state by one timestep dt using 4th-order Runge-Kutta integration
State rk4(State s, Params p, double dt) {

    State k1, k2, k3, k4;

    // k1 — derivatives at current state
    k1.theta1 = dtheta1(s, p);
    k1.omega1 = domega1(s, p);
    k1.theta2 = dtheta2(s, p);
    k1.omega2 = domega2(s, p);

    // k2 — derivatives at midpoint estimated using k1
    State s2;

    s2.theta1 = s.theta1 + 0.5 * dt * k1.theta1;
    s2.omega1 = s.omega1 + 0.5 * dt * k1.omega1;
    s2.theta2 = s.theta2 + 0.5 * dt * k1.theta2;
    s2.omega2 = s.omega2 + 0.5 * dt * k1.omega2;
    k2.theta1 = dtheta1(s2, p);
    k2.omega1 = domega1(s2, p);
    k2.theta2 = dtheta2(s2, p);
    k2.omega2 = domega2(s2, p);

    // k3 — derivatives at midpoint estimated using k2
    State s3;

    s3.theta1 = s.theta1 + 0.5 * dt * k2.theta1;
    s3.omega1 = s.omega1 + 0.5 * dt * k2.omega1;
    s3.theta2 = s.theta2 + 0.5 * dt * k2.theta2;
    s3.omega2 = s.omega2 + 0.5 * dt * k2.omega2;
    k3.theta1 = dtheta1(s3, p);
    k3.omega1 = domega1(s3, p);
    k3.theta2 = dtheta2(s3, p);
    k3.omega2 = domega2(s3, p);

    // k4 — derivatives at endpoint estimated using k3
    State s4;

    s4.theta1 = s.theta1 + dt * k3.theta1;
    s4.omega1 = s.omega1 + dt * k3.omega1;
    s4.theta2 = s.theta2 + dt * k3.theta2;
    s4.omega2 = s.omega2 + dt * k3.omega2;
    k4.theta1 = dtheta1(s4, p);
    k4.omega1 = domega1(s4, p);
    k4.theta2 = dtheta2(s4, p);
    k4.omega2 = domega2(s4, p);
    // final weighted combination
    State result;
    result.theta1 = s.theta1 + (dt/6.0) * (k1.theta1 + 2.0*k2.theta1 + 2.0*k3.theta1 + k4.theta1);
    result.omega1 = s.omega1 + (dt/6.0) * (k1.omega1 + 2.0*k2.omega1 + 2.0*k3.omega1 + k4.omega1);
    result.theta2 = s.theta2 + (dt/6.0) * (k1.theta2 + 2.0*k2.theta2 + 2.0*k3.theta2 + k4.theta2);
    result.omega2 = s.omega2 + (dt/6.0) * (k1.omega2 + 2.0*k2.omega2 + 2.0*k3.omega2 + k4.omega2);
    return result;
}
double total_energy(State s, Params p){
    return 0.5*p.m1*p.l1*p.l1*s.omega1*s.omega1 + 0.5*p.m2*(p.l1*p.l1*s.omega1*s.omega1 
        + p.l2*p.l2*s.omega2*s.omega2 + 2.0*p.l1*p.l2*s.omega1*s.omega2*std::cos(s.theta1 - s.theta2))
        -(p.m1 + p.m2)*p.g*p.l1*std::cos(s.theta1)- p.m2*p.g*p.l2*std::cos(s.theta2);
}
int main(){

    State start;
    start.theta1 = M_PI / 3.0;
    start.omega1 = 0.0;
    start.theta2 = M_PI / 3.0;
    start.omega2 = 0.0;

    Params parameter;
    parameter.l1 = 1.0;
    parameter.l2 = 1.0;
    parameter.m1 = 1.0;
    parameter.m2 = 1.0;
    parameter.b1 = 0.0;
    parameter.b2 = 0.0;
    parameter.g  = 9.81;

    int N=10000;

    for(int i=0;i<N;i++){
        start = rk4(start,parameter,1e-3);
        std::cout << start.theta1 << ", " << start.theta2 << ", " << total_energy(start, parameter) << std::endl;
    }

}






















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




double dtheta1(State s, Params p)  {
    return s.omega1;
}

double dtheta2(State s, Params p){
    return s.omega2;
}
double domega1(State s, Params p) {
    double dtheta = s.theta1 - s.theta2;
    double num = -p.m2 * p.g * std::sin(s.theta1 - 2.0*s.theta2)
                 - 2.0 * std::sin(dtheta) * (p.m2 * p.l2 * s.omega2*s.omega2 + p.l1 * s.omega1*s.omega1 * std::cos(dtheta))
                 - (p.m1 + p.m2) * p.g * std::sin(s.theta1);
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






























#include "parameters.h"
#include "vector.h"
#include <array>
#include <QVector>
#include <QDebug>

#include <cassert>
#include "matr_functions.h"

//nu1, nu2, nu3, x, y, theta

Vector<6> rightpart_custom(double t, Vector<6> x, Vector<3> control_minus, Vector<3> control_plus, double t_sw,
                           std::array<double, 3> alpha, std::array<double, 3> beta,
                           std::array<double, 3> delta, double Delta, double Lambda)
{
    Vector<6> rez;
    Vector<3> u;


    if (t < t_sw)
        u = control_minus;
    else
        u = control_plus;

    auto A_rev = matr_rev(A_matr(alpha, beta, delta, Delta, Lambda));
    auto Sigma = Sigma_matr(alpha, beta, delta, Delta, Lambda);


    auto second_matr = scal_mult(matr_mult(A_rev, matr_T(Sigma)), parameters::c1);
    auto third_matr = scal_mult(matr_mult(A_rev, matr_mult(matr_T(Sigma), Sigma)), -parameters::c2);

    auto rightpart_expr = vec_mult(A_rev, {parameters::symmetrical:: m * x[1] * x[2] / Lambda,
                                           - parameters::symmetrical:: m * x[0] * x[2] / Lambda, 0}) +
                            vec_mult(second_matr, u) + vec_mult(third_matr, {x[0], x[1], x[2]});

    rez[0] = rightpart_expr[0];
    rez[1] = rightpart_expr[1];
    rez[2] = rightpart_expr[2];
    rez[3] = x[0] * cos(x[5]) - x[1] * sin(x[5]);
    rez[4] = x[0] * sin(x[5]) + x[1] * cos(x[5]);
    rez[5] = x[2] / Lambda;

    return rez;
}

Vector<6> DOPRI8_custom(double t_left, double t_right, Vector<6> initial_values,
                    Vector<3> control_minus, Vector<3> control_plus, double t_sw,
                    std::array<double, 3> alpha, std::array<double, 3> beta,
                    std::array<double, 3> delta, double Delta, double Lambda)
{
    double h = (t_right - t_left) / 1e7;
    double h_new;
    bool switch_flag = false;
    bool last_flag = false;

    double tl = t_left;
    Vector<6> xl = initial_values;

    Vector<6> stepx;

    Vector<6> errx;
    double err, coef;
    double coefmax = 5;

    Vector<6> k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13;

    while (tl + h < t_right || last_flag)
    {

        //switch point
        if (tl < t_sw && tl + h > t_sw && !switch_flag)
        {
            h = t_sw - tl;
            switch_flag = true;
        }

        k1 = rightpart_custom(tl, xl, control_minus, control_plus, t_sw, alpha, beta, delta, Delta, Lambda);

        k2 = rightpart_custom(tl + h / 18, xl + h * k1 / 18, control_minus, control_plus, t_sw, alpha, beta, delta, Delta, Lambda);

        k3 = rightpart_custom(tl + h / 12, xl + h * (k1  / 48 + k2 / 16) , control_minus, control_plus, t_sw, alpha, beta, delta, Delta, Lambda);

        k4 = rightpart_custom(tl + h / 8, xl + h * (k1  / 32 + k3 * 3 / 32) , control_minus, control_plus, t_sw, alpha, beta, delta, Delta, Lambda);

        k5 = rightpart_custom(tl + h * 5 / 16, xl + h * (k1  * 5 / 16 + k3 * (-75) / 64 + k4 * 75 / 64) , control_minus, control_plus, t_sw, alpha, beta, delta, Delta, Lambda);

        k6 = rightpart_custom(tl + h * 3 / 8, xl + h * (k1  * 3 / 80 + k4 * 3 / 16 + k5 * 3 / 20) , control_minus, control_plus, t_sw, alpha, beta, delta, Delta, Lambda);

        k7 = rightpart_custom(tl + h * 59. / 400., xl + h * (k1  * (29443841. / 614563906.) + k4 * (77736538. / 692538347.) + k5 * ((-28693883.) / 1125000000.) +
                k6 * (23124283. / 1800000000.) ), control_minus, control_plus, t_sw, alpha, beta, delta, Delta, Lambda);

        k8 = rightpart_custom(tl + h * 93. / 200., xl + h * (k1  * (16016141. / 946692911.) + k4 * (61564180. / 158732637.) + k5 * (22789713. / 633445777.) +
                k6 * (545815736. / 2771057229.) + k7 * ((-180193667.) / 1043307555.) ), control_minus, control_plus, t_sw, alpha, beta, delta, Delta, Lambda);

        k9 = rightpart_custom(tl + h * (5490023248. / 9719169821.), xl + h * (k1  * (39632708. / 573591083.) + k4 * ((-433636366.) / 683701615.) + k5 * ((-421739975.) / 2616292301.) +
                k6 * (100302831. / 723423059.) + k7 * (790204164. / 839813087.) + k8 * (800635310. / 3783071287.) ), control_minus, control_plus, t_sw, alpha, beta, delta, Delta, Lambda);

        k10 = rightpart_custom(tl + h * 13 / 20, xl + h * (k1  * (246121993. / 1340847787.) + k4 * ((-37695042795.) / 15268766246.) + k5 * ((-309121744.) / 1061227803.) +
                k6 * ((-12992083.) / 490766935.) + k7 * (6005943493. / 2108947869.) + k8 * (393006217. / 1396673457) + k9 * (123872331. / 1001029789.) ), control_minus, control_plus, t_sw, alpha, beta, delta, Delta, Lambda);

        k11 = rightpart_custom(tl + h * (1201146811. / 1299019798.), xl + h * (k1  * ((-1028468189.) / 846180014.) + k4 * (8478235783. / 508512852.) + k5 * (1311729495. / 1432422823.) +
                k6 * ((-10304129995.) / 1701304382.) + k7 * ((-48777925059.) / 3047939560.) + k8 * (15336726248. / 1032824649.) + k9 * ((-45442868181.) / 3398467696.) +
                k10 * (3065993473. / 597172653.) ), control_minus, control_plus, t_sw, alpha, beta, delta, Delta, Lambda);

        k12 = rightpart_custom(tl + h, xl + h * (k1  * (185892177. / 718116043.) + k4 * ((-3185094517.) / 667107341.) + k5 * ((-477755414.) / 1098053517.) +
                 k6 * ((-703635378.) / 230739211.) + k7 * (5731566787. / 1027545527.) + k8 * (5232866602. / 850066563.) + k9 * ((-4093664535.) / 808688257.) +
                k10 * (3962137247. / 1805957418.) + k11 * (65686358. / 487910083.) ), control_minus, control_plus, t_sw, alpha, beta, delta, Delta, Lambda);


        k13 = rightpart_custom(tl + h, xl + h * (k1  * (403863854. / 491063109.) + k4 * ((-5068492393.) / 434740067.) + k5 * ((-411421997.) / 543043805.) +
                 k6 * (652783627. / 914296604.) + k7 * (11173962825. / 925320556.) + k8 * ((-13158990841.) / 6184727034.) + k9 * (3936647629. / 1978049680.) +
                k10 * ((-160528059.) / 685178525.) + k11 * (248638103. / 1413531060.) ), control_minus, control_plus, t_sw, alpha, beta, delta, Delta, Lambda);


        stepx = h * (k1 * (14005451. / 335480064.) + k6 * ((-59238493.) / 1068277825.) + k7 * (181606767. / 758867731.) + k8 * (561292985. / 797845732.) +
                k9 * ((-1041891430.) / 1371343529.) + k10 * (760417239. / 1151165299.) + k11 * (118820643. / 751138087.) + k12 * ((-528747749.) / 2220607170.) + k13 / 4 );


        errx = h * (k1 * (13451932. / 455176623.) + k6 * ((-808719846.) / 976000145.) + k7 * (1757004468. / 5645159321.) + k8 * (656045339. / 265891186.) +
                k9 * ((-3867574721.) / 1518517206.) + k10 * (465885868. / 322736535.) + k11 * (53011238. / 667516719.) + k12 * 2 / 45 );


        err = (stepx - errx).euclidDistance();

        //step correction section
        if (err < 1e-15)
        {
            tl += h;
            xl += stepx;

            coefmax = 5;

            h *= 2;
        }
        else
        {
            coef = pow(precision::DOPRI8_error_EPS / err, 1.0 / (8 + 1));

            if (coef > coefmax)
                coef = coefmax;
            else if (coef < 0.2)
                coef = 0.2;

            h_new = 0.9 * h * coef;


            if (precision::DOPRI8_error_EPS > err)
            {
                tl += h;
                xl += stepx;

                coefmax = 5;
            }
            else
                coefmax = 1;

            h = h_new;
        }

        //last step correction
        if (tl + h >= t_right && last_flag == false)
        {
            h = t_right - tl;
            last_flag = true;
        }
        else
            last_flag = false;
    }

    return xl;
}

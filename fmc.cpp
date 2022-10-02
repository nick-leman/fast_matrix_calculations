#include <iostream>
#include <conio.h>
#include <math.h>
#define SIZE 20 //SIZE = ðàíã ìàòðèö³ äëÿ ïðîì³æíèõ îá÷èñëåíü; ðàíã ìàòðèö³ = SIZE/2-1 (ìîæíà ñòàâèòè á³ëüøå í³æ ïîòð³áíî)
using namespace std;
class matrix {   
private:
    double* p_m;
    int n;
public:
    matrix(int count_n);//default and size
    matrix(const matrix& m);//values constructor
    matrix operator+(const matrix my);
    matrix operator-(const matrix my);
    matrix operator*(const matrix my);
    matrix operator/(const matrix my);
    void   operator=(const matrix my);
    matrix transpon(void) {
        double a[SIZE][SIZE] = { 0 };
        matrix m_ret(n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                a[j][i] = p_m[i + j * n];
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                m_ret.p_m[i + j * n] = a[i][j];
        return m_ret;
    }
    matrix inverse(void) {

        double d;
        int i, j, k;
        double a[SIZE][SIZE] = { 0 };
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= n; j++) {
                a[i][j] = p_m[(j - 1) + (i - 1) * n];
            }



        for (i = 1; i <= n; i++)
            for (j = 1; j <= 2 * n; j++)
                if (j == (i + n))
                    a[i][j] = 1;

        /************** partial pivoting **************/
        for (i = n; i > 1; i--)
        {
            if (a[i - 1][1] < a[i][1])
                for (j = 1; j <= n * 2; j++)
                {
                    d = a[i][j];
                    a[i][j] = a[i - 1][j];
                    a[i - 1][j] = d;
                }
        }
        /********** reducing to diagonal  matrix ***********/

        for (i = 1; i <= n; i++)
        {
            for (j = 1; j <= n * 2; j++)
                if (j != i)
                {
                    d = a[j][i] / a[i][i];
                    for (k = 1; k <= n * 2; k++)
                        a[j][k] -= a[i][k] * d;
                }
        }
        /************** reducing to unit matrix *************/
        for (i = 1; i <= n; i++)
        {
            d = a[i][i];
            for (j = 1; j <= n * 2; j++)
                a[i][j] = a[i][j] / d;
        }



        matrix m_ret(n);
        for (i = 1; i <= n; i++)
        {
            for (j = n + 1; j <= n * 2; j++) {
                if (isnan(a[i][j])||isinf(a[i][j]))  throw invalid_argument("Invalid operation: inverse matrix doesn`t exist.");
                m_ret.p_m[(i - 1) + (j - n - 1) * n] = a[i][j];
            }
        }

        return m_ret.transpon();
    }
    double det(void) {

        double det = 0;
        int per = 0;
        double d;
        int i, j, k;
        double a[SIZE][SIZE] = { 0 };
        for (int i = 1; i <= n; i++)
            for (int j = 1; j <= n; j++) {
                a[i][j] = p_m[(j - 1) + (i - 1) * n];
            }



        /************** partial pivoting **************/
        for (i = n; i > 1; i--)
        {
            if (a[i - 1][1] < a[i][1]) {
                for (j = 1; j <= n * 2; j++)
                {
                    d = a[i][j];
                    a[i][j] = a[i - 1][j];
                    a[i - 1][j] = d;
                }
                per++;
            }

        }

        /********** reducing to diagonal  matrix ***********/
        for (i = 1; i <= n; i++)
        {
            for (j = 1; j <= n * 2; j++)
                if (j != i)
                {
                    d = a[j][i] / a[i][i];
                    for (k = 1; k <= n * 2; k++)
                        a[j][k] -= a[i][k] * d;
                }
        }
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {

                if (i == j) det = (det == 0) ? a[i][j] : det * a[i][j];
            }
        det = (isnan(det)) ? 0 : det;
        return det * (pow(-1, per));
    }
    void   In_put();
    void   Out_put();
    ~matrix();
};
matrix  matrix::operator*(const matrix my)
{
    matrix m_ret(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double sum = 0;
            for (int k = 0; k < n; k++) sum += p_m[i * n + k] * my.p_m[k * n + j];
            m_ret.p_m[i * n + j] = sum;
        }
    }
    return m_ret;
};
matrix  matrix::operator-(const matrix my)
{
    int my_n = (n * n);
    matrix m_ret(n);
    for (int i = 0; i < my_n; i++)
        m_ret.p_m[i] = p_m[i] - my.p_m[i];
    return m_ret;
};
matrix matrix::operator+(const matrix my)
{
    int my_n = (n * n);
    matrix m_ret(n);
    for (int i = 0; i < my_n; i++)
        m_ret.p_m[i] = p_m[i] + my.p_m[i];
    return m_ret;
};
void matrix::operator=(const matrix my)
{
    int my_n = (n * n);
    for (int i = 0; i < my_n; i++)
        p_m[i] = my.p_m[i];
};
matrix matrix::operator/(const matrix my)
{
    int my_n = (n * n);
    matrix m_ret1(n);
    matrix m_ret2(n);
    for (int i = 0; i < my_n; i++) {
        m_ret1.p_m[i] = p_m[i];

    };
    for (int i = 0; i < my_n; i++) {
        m_ret2.p_m[i] = my.p_m[i];


    };
    m_ret2 = m_ret2.inverse();

    m_ret1 = m_ret1 * m_ret2;
    return  m_ret1;
};
matrix::matrix(int count_n)
{
    n = count_n;
    int my_n = (n * n);
    p_m = new double[my_n];
    for (int i = 0; i < my_n; i++)
        p_m[i] = 0;
};
matrix::matrix(const matrix& m) {
    n = m.n;
    int my_n = (n * n);
    p_m = new double[my_n];
    for (int i = 0; i < my_n; i++)
        p_m[i] = m.p_m[i];
};
matrix::~matrix()
{
    delete p_m;
};
void matrix::Out_put()
{
    int y = 0, out_ch = 0;
    for (int i = 0; i < n; i++)
    {
        for (int k = 0; k < y; k++)
        {
            cout.width(7);
        }
        for (int j = 0; j < n; j++)
        {
            cout.width(10);
            cout << round(p_m[out_ch] * 1000000) / 1000000;
            out_ch++;
        }
        cout << "\n";
        y++;
    }
};
void matrix::In_put()
{
    int y = 0, in_ch = 0;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            cout << "[" << i << "," << j << "]: ";
            cin >> p_m[in_ch];
            in_ch++;
        }
        y++;
    }
};
int main()
{
    int size;
    cout << "Set square matrix size: ";
    cin >> size;
    matrix mat1(size), mat2(size), mat(size);
    cout << "\nMatrix 1: \n";
    mat1.In_put();
    cout << "\nMatrix 2: \n";
    mat2.In_put();

    cout << "\n Matrix 1 + Matrix 2 : \n";
    mat = mat1 + mat2;
    mat.Out_put();
    cout << "\n";
    cout << "\n Matrix 1 - Matrix 2 : \n";
    mat = mat1 - mat2;
   
    mat.Out_put();
    cout << "\n";
    cout << "\n  Matrix 1 * Matrix 2  :\n";
    mat = mat1 * mat2;
    mat.Out_put();
    cout << "\n";
    cout << "\n Matrix 1 determ: ";
    cout << mat1.det() << "\n";
    cout << "\n";
    
    cout << "\n Matrix 2 inversed:\n";
    mat = mat2.inverse();
    mat.Out_put();
    cout << "\n";
    cout << "\n Matrix 1 transposed:\n";
    mat = mat1.transpon();
    mat.Out_put();
    cout << "\n Matrix 1 / Matrix 2 :\n";
    mat = mat1 / mat2;
    mat.Out_put();
    _getch();
};

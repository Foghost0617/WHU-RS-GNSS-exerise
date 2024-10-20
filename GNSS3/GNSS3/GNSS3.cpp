// GNSS3.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include<math.h>
#include <iomanip>
#include <cmath>
using namespace std;
extern double miu = 3.986004415e14;
extern double radv = 7.292115147e-5;
struct satellite {
    string name;//1 0
    int year,month,day,hour,min,sec;//1 1
    double af0;//1 2
    double af1;//1 3
    double af2;//1 4
    double temp21;//2 1
    double Crs;//2 2
    double deltan;//2 3
    double M0;//2 4
    double Cuc;//3 1
    double e;//3 2
    double Cus;//3 3
    double sqrtA;//3 4
    double toe;//4 1
    double Cic;//4 2
    double omega0;//4 3
    double Cis;//4 4
    double i0;//5 1
    double Crc;//5 2
    double omega;//5 3
    double omegaDot;//5 4
    double IDOT;//6 1
    double temp62, temp63, temp64;
    double temp71, temp72;
    double TGD;//7 3

    void print() const {
        std::cout <<name<< " "<<year<<" " <<month<<" "<<day<<" "<<hour<<" "<<min<<" "<<sec<<" "<<af0<<" "<<af1<<" "<<temp21 <<" "<< Crs << std::endl;
    }
};
struct SatTime {
    string c;
    int year, month, day, hour, min, sec;
};
struct SatTimeXYZ {
    string name;
    int year, month, day, hour, min, sec;
    double X, Y, Z;
};
struct user {
    string name;
    int year, month, day, hour, min, sec;
};

void Read1(const char* pFileName, vector< satellite>& Data) {
    string dataline;
    std::ifstream in(pFileName);
    if (!in.is_open()) {
        cout << "Failed" << endl;
    }
    std::string line;
    int skip = 6;
    for (int i = 0; i < skip; ++i) {
        getline(in, line);
    }//跳文件头    
    //vector<satellite1> dataList; //存多个结构体
    while (getline(in, line)) {
        if (line.empty() || line[0] == 'R') {
            for (int i = 0; i < 3; ++i) {
                getline(in, line);
            }
            continue;
        }
        else if (line.empty() || line[0] == 'C') {
            for (int i = 0; i < 7; ++i) {
                getline(in, line);
            }
            continue;
        }
        else if (line.empty() || line[0] == 'G') {
            for (int i = 0; i < 7; ++i) {
                dataline = dataline + line;
                //if (i == 7) { break; }
                getline(in, line);
            }
            //cout << dataline << endl;
            istringstream iss(dataline);//按空格回车分割
            satellite g;
            iss >> g.name >> g.year >> g.month >> g.day >> g.hour >> g.min >> g.sec >> g.af0 >> g.af1 >> g.af2 >> g.temp21 >> g.Crs >> g.deltan >> g.M0 >> g.Cuc >> g.e >> g.Cus >> g.sqrtA >> g.toe >> g.Cic >> g.omega0 >> g.Cis >> g.i0 >> g.Crc >> g.omega >> g.omegaDot >> g.IDOT;
            Data.push_back(g);
            dataline = " ";
        }
    }
    in.close();
}
void Read2(const char* pFileName,vector<SatTimeXYZ>&SatPos) {
    //string dataline;
    std::ifstream in(pFileName);
    if (!in.is_open()) {
        cout << "Failed" << endl;
    }
    string line;
    int skip = 22;
    for (int i = 0; i < skip; ++i) {
        getline(in, line);
    }//跳文件头    
    while (getline(in, line)) {
        if (line[0] == '*') {
            istringstream iss(line);
            SatTime time;
            iss >> time.c >> time.year >> time.month >> time.day >> time.hour >> time.min >> time.sec;
            getline(in, line);
            if (line[0] == 'P' && line[1] == 'G') {
                for(int i=0;i<=31;++i){
                    istringstream iss(line);
                    SatTimeXYZ g;
                    iss >> g.name >> g.X >> g.Y >> g.Z;
                    g.name=g.name.erase(0, 1);
                    g.year = time.year; g.month = time.month; g.day = time.day; g.hour = time.hour; g.min = time.min; g.sec = time.sec;
                    SatPos.push_back(g);
                    if (i == 31) 
                    { break; }
                    else 
                    { getline(in, line); }
                } 
            }
        }
       
    }
}
void Compare(user& UserData, vector< satellite>& Data1, vector<SatTimeXYZ>& Data2) {
    int sign = 0;
    //t-toe
    int delta_t=0;
    for (size_t i = 0; i < Data1.size(); ++i) {
        if (UserData.name == Data1[i].name && UserData.year == Data1[i].year && UserData.month == Data1[i].month && UserData.day == Data1[i].day && UserData.hour >= Data1[i].hour) {//判断卫星号年月日相同
            int sec1 = UserData.hour * 3600 + UserData.min * 60 + UserData.sec;
            int sec2 = Data1[i].hour * 3600 + Data1[i].min * 60 + Data1[i].sec;
            if (sec1 >= sec2) { 
                delta_t = sec1 - sec2;//算秒差
                sign = i;//标记
                break;//找到匹配后退循环
            }
        }
    }
    //cout << "delta_t " << delta_t << endl;

    //M
    double M = 0.0;
    double n0 = sqrt(miu) / (Data1[sign].sqrtA * Data1[sign].sqrtA * Data1[sign].sqrtA);
    double n = n0 + Data1[sign].deltan;
    M = Data1[sign].M0 + n * delta_t;
    //cout << M << endl;

    //偏近点角E=M+e*sinE    
    double E0 = M;
    double E = 0.0;
    double delta=0;
    do {
        E = M + Data1[sign].e * sin(E0);
        delta = E - E0;
        E0 = E;
    } while (abs(delta) >= 1e-12);
    //cout << setiosflags(ios::scientific) << setprecision(16) <<"E "<< E << endl;
    //cout << "E: " << E << endl;

    //真近点角
    double tempY = sqrt(1 - Data1[sign].e * Data1[sign].e) * sin(E);
    double tempX = cos(E) - Data1[sign].e;
    double f = atan2(tempY, tempX);
    
    //cout <<"f " << f << endl;
    
    //升交角距/未改正
    double u1 = Data1[sign].omega + f;
    //cout << u1 << endl;

    //卫星向径/未改正
    double r1 = (Data1[sign].sqrtA) * (Data1[sign].sqrtA) * (1 - Data1[sign].e * cos(E));
    //cout << r1 << endl;

    //摄动改正项
    double u_ =Data1[sign].Cuc * cos(2 * u1) + Data1[sign].Cus * sin(2 * u1);
    double r_ = Data1[sign].Crc * cos(2 * u1) + Data1[sign].Crs * sin(2 * u1);
    double i_ = Data1[sign].Cic*cos(2*u1)+ Data1[sign].Cis * sin(2 * u1);
    //cout << "u_ " << u_ << " r_ " << r_ << " i_ " << i_ << endl;

    //摄动改正
    double u = u1 + u_;
    double r = r1 + r_;
    double i = Data1[sign].i0 + Data1[sign].IDOT * delta_t + i_;
    //cout << "u " << u << " r " << r << " i " << i << endl;

    //计算卫星在轨道平面坐标系中的位置
    double x = r * cos(u);
    double y = r * sin(u);
    //cout << x << " " << y << endl;

    //升交点经度L
    double L = Data1[sign].omega0 + Data1[sign].omegaDot * delta_t - radv*(delta_t + Data1[sign].toe);
    //cout <<"L " << L << endl;

    //计算瞬时地球坐标系下坐标
    double X = x * cos(L) - y * cos(i) * sin(L);
    double Y = x * sin(L) + y * cos(i) * cos(L);
    double Z = y * sin(i);
    cout << "卫星" << UserData.name << "在" << UserData.year << " " << UserData.month << " " << UserData.day << " " << UserData.hour << " " << UserData.min << " " << UserData.sec << "的坐标" << endl;
    cout << fixed << setprecision(9)<< "X：" << X << " Y：" << Y << " Z：" << Z <<endl;
    
    ////协议 
    //double X = X1 + Xp * Z1;
    //double Y = Y1 - Yp * Z1;
    //double Z = Z1+ Yp*Y1-X1*Xp;
    //cout <<"X: " << X << " Y: " << Y << " Z: " << Z << endl;

    double delta_X=0, delta_Y=0, delta_Z = 0;
    for (size_t i = 0; i < Data2.size(); ++i) {
        if (UserData.name == Data2[i].name && UserData.year == Data2[i].year && UserData.month == Data2[i].month && UserData.day == Data2[i].day && UserData.hour == Data2[i].hour && UserData.min == Data2[i].min && UserData.sec == Data2[i].sec) {
            delta_X = X - 1000 * Data2[i].X;
            delta_Y = Y - 1000 * Data2[i].Y;
            delta_Z = Z - 1000 * Data2[i].Z;
            break;
        }
    }
    cout << "与精密星历较差：" << endl;
    cout << "delta_X: " << delta_X << " delta_Y: " << delta_Y << " delta_Z: " << delta_Z << endl;
}

int main()
{
    const char* pFileName1 = "SJ111708S.23P";
    vector<satellite> Data1;
    Read1(pFileName1, Data1);
    const char* pFileName2 = "igv22671_18.sp3";
    vector<SatTimeXYZ>Data2;
    Read2(pFileName2, Data2);

    string UserInput1="G17 2023 6 19 2 00 0";
    string UserInput2 = "G17 2023 6 19 2 30 0";
    string UserInput3 = "G17 2023 6 19 3 00 0";
    
    user input1;
    istringstream iss1(UserInput1);
    iss1 >> input1.name >> input1.year >> input1.month >>input1.day>> input1.hour >> input1.min >> input1.sec;
    Compare(input1, Data1, Data2);
    cout << endl;
    user input2;
    istringstream iss2(UserInput2);
    iss2 >> input2.name >> input2.year >> input2.month >> input2.day >> input2.hour >> input2.min >> input2.sec;
    Compare(input2, Data1, Data2);
    cout << endl;
    user input3;
    istringstream iss3(UserInput3);
    iss3 >> input3.name >> input3.year >> input3.month >> input3.day >> input3.hour >> input3.min >> input3.sec;
    Compare(input3, Data1, Data2);
    cout << endl;

}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件

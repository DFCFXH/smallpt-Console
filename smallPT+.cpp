/*路径追踪的基本思想是从视点发出一条光线,光线与物体表面相交时根据表面的材质属性继续采样一个方向,发出另一条光线,如此
迭代,直到光线打到光源上(或逃逸出场景),然后用蒙特卡洛的方法,计算其贡献,作为像素的颜色值.
路径追踪会避开贡献小的路径,而在贡献大的路径附近做更多局部的探索.
路径追踪=光线追踪+蒙特卡洛方法*/
#define _CRT_SECURE_NO_WARNINGS//禁用fopen的警告
//#define _OPENMP
#include <math.h>//smallpt,一个小型路径追踪器
#include <stdlib.h> 
#include <stdio.h>
#include <omp.h>
#include <iostream>
#include <format>
#include <conio.h>
#include <sstream>
#include <Windows.h>
#include <cmath>
using namespace std;
//#define double float//使用单精度浮点数以节省内存
#define float2 double
#define M_PI 3.1415926525//使用常量M_PI表示圆周率
//#define IMAGE_PATHNAME "image.ppm"
double erand48(unsigned short xsubi[3]) {//随机数
    return (double)rand() / (double)RAND_MAX;
}
struct Vec {//三维向量
    double x, y, z;//位置或颜色都能使用
    Vec(double x_ = 0, double y_ = 0, double z_ = 0) { x = x_; y = y_; z = z_; }//构造函数,x,y,z默认为零
    Vec operator+(const Vec& b) const { return Vec(x + b.x, y + b.y, z + b.z); }
    Vec operator-(const Vec& b) const { return Vec(x - b.x, y - b.y, z - b.z); }
    Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }//向量乘标量
    Vec mult(const Vec& b) const { return Vec(x * b.x, y * b.y, z * b.z); }//向量乘向量,用于颜色计算
    Vec& norm() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }//归一化
    double dot(const Vec& b) const { return x * b.x + y * b.y + z * b.z; }//点乘
    Vec operator%(Vec& b) { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); }//叉乘
    Vec operator/(const Vec& b) const { return Vec(x / b.x, y / b.y, z / b.z); }
    double DIV(const double& b) const { return x /b + y / b + z / b; }//向量除以标量
};
Vec VecRotate(Vec vec, double angle, Vec axis) {
    double c = cos(angle);
    double s = sin(angle);
    double t = 1 - c;
    double x = vec.x, y = vec.y, z = vec.z;
    double newX = (t * axis.x * axis.x + c) * x + (t * axis.x * axis.y - s * axis.z) * y + (t * axis.x * axis.z + s * axis.y) * z;
    double newY = (t * axis.x * axis.y + s * axis.z) * x + (t * axis.y * axis.y + c) * y + (t * axis.y * axis.z - s * axis.x) * z;
    double newZ = (t * axis.x * axis.z - s * axis.y) * x + (t * axis.y * axis.z + s * axis.x) * y + (t * axis.z * axis.z + c) * z;
    return Vec(newX, newY, newZ);
}
struct Ray { Vec o, d; Ray(Vec o_, Vec d_) : o(o_), d(d_) {} };//光线结构体
enum obj_type { S, P, T,BVH };//材质类型
struct Obj {
    double rad, ref, diff, spec, refr, refr_nt,w,h,size;//半径,反射,漫反射,镜面反射,折射概率,折射率,宽,高
    Vec p, e, c, i, n,p1,p2,p3;//位置,自发光,颜色(这里的颜色并不是0-255,而是0-1),杂质,法线
    obj_type type;// 物体类型
    Obj(double rad_, double ref_, double diff_, double spec_, double refr_, double refr_nt_,Vec i_,Vec p_, Vec e_, Vec c_,
        obj_type refl_,double w_,double h_,Vec n_,double size_,Vec p1_,Vec p2_,Vec p3_) ://构造函数
        rad(rad_), ref(ref_), diff(diff_), spec(spec_), refr(refr_), refr_nt(refr_nt_),i(i_), p(p_), e(e_), c(c_), type(refl_),w(w_),h(h_),n(n_)
    ,size(size_),p1(p1_), p2(p2_), p3(p3_) {}
    double intersect_sphere(const Ray& r) const { //计算射线原点与球体之间的交点的距离,如果没有交点返回0
        Vec op = p - r.o; //光源指向球心的一条向量
        double t, eps = 1e-4, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
        //eps是一个很小的量,代指0
        //t是射线与交点之间的距离,是方程t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0的解"."表示点乘
        //b是光源指向球心的向量和光源方向向量的夹角的余弦值
        //det没有解(<0)则没有相交^^^|op-t|=|r|         (op-t)^2-r^2=0
        if (det < 0) return 0; else det = sqrt(det);
        return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
        //返回t值,如果为0,则表示不相交
        //选择其中较小且大于0的解,如果det<0则无解,表示不相交,返回0
    };
    double intersect_plane(const Ray& r) const {
        double t = n.dot(p - r.o) / n.dot(r.d);
        double dn = n.dot(r.d);

        // 如果光线平行于或位于平面后面，则没有相交
        if (fabs(dn) < 1e-6 || dn >= 0) {
            return 0;
        }

        // 交点
        Vec intersection = r.o + r.d * t;

        // 检查交点是否在平面范围内
        double half_width = w / 2.0;
        double half_height = h / 2.0;
        double z_distance = intersection.z - p.z; // 计算交点到平面中心的Z轴距离

        if (intersection.x < p.x - half_width || intersection.x > p.x + half_width ||
            intersection.y < p.y - half_height || intersection.y > p.y + half_height ||
            z_distance < 0 || z_distance > h) { // 检查交点在Z轴方向上是否在范围内
            return 0;
        }

        // 返回相交点与光线原点之间的距离
        return t;
    };
    double intersect_triangle(const Ray& r) const { //计算射线原点与球体之间的交点的距离,如果没有交点返回0
        /*Vec p2_ = p2 - p1;
        Vec p3_ = p3 - p1;
        Vec p1 = p;
        Vec p2 = p1 + p2_;
        Vec p3 = p1 + p3_;
        Vec nl = n;
        if (n.dot(r.d) > 0) nl=n * -1;
        if (fabs(nl.dot(r.d)) < -1e4) return 0;
        double t = (nl.dot(p1) - r.o.dot(nl)) / r.d.dot(nl);
        if (t < 0.0005f) return 0;
        Vec P = r.o +r.d*t;
        Vec a, b;
        a = p2 - p1;
        b = P - p1;
        Vec c1 = a%b;
        a = p3 - p2;
        b = P - p2;
        Vec c2 = a % b;
        a = p1 - p3;
        b = P - p3;
        Vec c3 = a % b;
        if (c1.dot(n) < 0 || c2.dot(n) < 0 || c3.dot(n) < 0) return 0;
        return t;*/
        const double EPSILON = 0.0000001; //要比较的小值
        Vec p1_ = p1 * size + p;
        Vec p2_ = p2 * size + p;
        Vec p3_ = p3 * size + p;
        Vec p2 = p + p2_;
        Vec p3 = p + p3_;
        Vec p1 = p + p1_;
        Vec edge1 = p2 - p1;
        Vec edge2 = p3 - p1;
        Vec rd = r.d;
        Vec h = rd % edge2;
        double a = edge1.dot(h);
        Vec AB = p2 - p1, AC = p3 - p1, n0 = AB % AC;
        Vec n = n0.norm();
        double tdn = n.dot(r.d);
        if (a > -EPSILON && a < EPSILON || tdn >= 0)
            return 0.0; //射线平行于三角形

        double f = 1.0 / a;
        Vec s = r.o - p1;
        double u = f * s.dot(h);

        if (u < 0.0 || u > 1.0)
            return 0.0; //交点在三角形之外

        Vec q = s % edge1;
        double v = f * r.d.dot(q);

        if (v < 0.0 || u + v > 1.0)
            return 0.0; //交点在三角形之外

        double t = f * edge2.dot(q);

        if (t > EPSILON)
            return t; //找到交点

        return 0.0; //未找到交点
    };
};
/*半径,反射概率,漫反射概率,镜面反射概率,折射概率,折射率,杂质,位置,自发光,颜色,类型,
宽,高,法线,大小,三点位置*/
Obj scenes[] = {
      Obj(13,100,100,0,0,1.5,Vec(),Vec(),Vec(),Vec(1,.45,1) * .999, P,
  1000,1000,Vec(0,1,0)
  ,1,Vec(1,1,1),Vec(2,2,2),Vec(3,3,3)),//
    Obj(600,100,100,0,0,1.5,Vec(),Vec(0,800,0),Vec(24,24,24),Vec(1,.45,1) * .999, S,
  1000,1000,Vec(0,1,0)
  ,1,Vec(1,1,1),Vec(2,2,2),Vec(3,3,3)),//
   Obj(30,100,100,0,0,1.5,Vec(),Vec(50,80,0),Vec(),Vec(1,.45,0) * .999, S,
  1000,1000,Vec(0,1,0)
  ,1,Vec(1,1,1),Vec(2,2,2),Vec(3,3,3)),//
   Obj(60,100,0,100,0,1.5,Vec(),Vec(170,80,0),Vec(),Vec(0,0,1) * .999, S,
  1000,1000,Vec(0,1,0)
  ,1,Vec(1,1,1),Vec(2,2,2),Vec(3,3,3)),//
   Obj(60,100,0,0,100,3,Vec(),Vec(300,80,0),Vec(),Vec(1,1,1) * .999, S,
  1000,1000,Vec(0,1,0)
  ,1,Vec(1,1,1),Vec(2,2,2),Vec(3,3,3)),//
};
//这两个函数是亮度和颜色计算的辅助函数
inline double clamp(double x) { return x < 0 ? 0 : x>1 ? 1 : x; }//对于经过递归叠加后大于1的设为1,小于0的设为0
inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }/*把0 - 1转换为rgb中的0 - 255, 设置了各1 / 2.2的调整值,
能让画面更亮*/
//光线求交
inline bool intersect(const Ray& r, double& t, int& id) {//t表示距离
    double n = sizeof(scenes) / sizeof(Obj), d, inf = t = 1e20;
    for (int i = int(n); i--;) {//求出最近交点
        if (scenes[i].type == S) {
            if ((d = scenes[i].intersect_sphere(r)) && d < t) {
                t = d;
                id = i;
            }
        }
        else if (scenes[i].type == P) {
            if ((d = scenes[i].intersect_plane(r)) && d < t) {
                t = d;
                id = i;
            }
        }
        else if (scenes[i].type == T) {
            if ((d = scenes[i].intersect_triangle(r)) && d < t) {
                t = d;
                id = i;
            }
        }
    }
    return t < inf;//检查t是否过大,如果t过大则没有相交,返回false,否则返回true
}
extern double ncg = 1;//空气折射率
int numSpheres = sizeof(scenes) / sizeof(Obj);
Vec radiance(const Ray& r, int depth, unsigned short* Xi,int E=1) {
    double t;// 相交距离
    int id = 0;// 相交对象的ID
    if (!intersect(r, t, id)) return Vec(); // 未相交则返回黑色
    const Obj& obj = scenes[id];//被击中的对象
    Vec x = r.o + r.d * t,  f = obj.c;
    Vec n,pn=obj.n;
    Vec AB = obj.p2 - obj.p1, AC = obj.p3 - obj.p1, n0_t = AB % AC;
    if (obj.type == S) {
        n = (x - obj.p).norm();
    }
    else if (obj.type == P) {
        n = (pn).norm();
    }
    else if (obj.type == T) {
        n = n0_t.norm();
    }
    Vec nl = n.dot(r.d) < 0 ? n : n * -1;
    int fanshe = rand() % 100;
    int fanshe2 = rand() % 100;
    int zheshe = rand() % 100;
    /*x为交点,n为球体法向量,nl用于修正法向量,如果球体法向量和光线方向向量的点积小于零,则法线变为相反方向
    (此时光线从内部发出),(从外部来的光线,法线向外;从内部来的光线,法线向内),f为球体颜色*/
    if (depth > 6) return Vec();//当递归深度大于6,返回黑色
    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y : f.z; /*获取RGB三个值里的最高值*/
    if (++depth > 5) {//递归达到5时有机会返回
        if (erand48(Xi) < p)
            f = f * (1 / p);
        else
            return obj.e;
    }
    if (fanshe <= obj.ref) {
        if (obj.diff >= fanshe2) {// 漫反射(在半球当中随即找一个方向,然后进行递归)
            double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
            //r1为随机选取的角度,范围是 0 到 2π 之间,r2是随机选择了一个距离(0-1),r2s是距离开方的结果
            Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1, 0) : Vec(1, 0, 0)) % w).norm(), v = w % u;
            //fabs()求浮点数绝对值
            /*根据法线构造了一个正交基, w与法线同向(在n维空间当中,由n个互相正交(垂直)的向量组成一个正交基)
            在w.x的绝对值>0.1的时候,u垂直于(0,1,0)和w的单位向量,否则是垂直于(1,0,0)和w的单位向量
            这样做的目的是当w.x等于或接近0时,可能会出现线性相关的情况*/
            Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();//反射光方向向量
            return obj.e + f.mult(radiance(Ray(x, d), depth, Xi));
        }/*全局光照方程,使用蒙特卡罗方法求解*/
        /*镜面反射,
         公式为 入射光线方向向量-2*法向量*入射光方向向量dot法向量=反射光线方向向量*/
        else if (obj.spec >= fanshe2) {
            Vec spec_t= r.d - n * 2 * n.dot(r.d);
            if (obj.type == S) {
                return obj.e + f.mult(radiance(Ray(x, spec_t), depth, Xi));
            }
            else if (obj.type==P){
                return obj.e + f.mult(radiance(Ray(x, spec_t), depth, Xi));
            }
        }
        //以下为折射
        Ray reflRay(x, r.d - n * 2 * n.dot(r.d));// 反射光线
        bool into = n.dot(nl) > 0; /* 入射的光线是否从外面进来?如果n和nl同向则光线从外边进来,否则光线从内部发出*/
        double nc = ncg, nt = obj.refr_nt, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;/*nc为空气的折射率,nt为玻璃的折射
        率,nnt是原介质和目标介质的比值,ddn是光线方向向量和法线nl的夹角余弦值,cos2t是cos(t)^2*/
        double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
        Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1, 0) : Vec(1, 0, 0)) % w).norm(), v = w % u;
        Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
        if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0)    /*全反射,当光线从较高 折射率 的 介质 进入到较低折射率的介质
            时,如果入射角大于某一临界角θc(光线远离 法线 )时,折射角将变得足够大,折射后的光线将不会离开介质*/
            return obj.e + f.mult(radiance(reflRay, depth, Xi));
        Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();//折射光线的方向
        return obj.refr >= fanshe2 ? obj.e + f.mult(radiance(Ray(x, tdir), depth, Xi))+obj.i  :zheshe<=50?
            obj.e+radiance(reflRay, depth, Xi)* (obj.spec)*0.9 : obj.e + f.mult(radiance(Ray(x, d), depth, Xi)) * (obj.diff)*0.9;
    }
    else {
        return Vec();
    }
}
#define jz 0.01
Vec DeNoisy(int x, int y, int i, int w, int h, int samps, Vec c[], int i1=0) {
    Vec yansehuancun;
    double p = c[i].x > c[i].y && c[i].x > c[i].z ? c[i].x : c[i].y > c[i].z ? c[i].y : c[i].z;
    double a = 0;
    for (int j = 1; j <= i1; j++) {
        if (c[i].x <= jz && c[i].y <= jz && c[i].z <= jz && x > 0 && y > 0 && samps <= 100) {//降噪
            if (y > 0 + j && i - w - j <= w * h && i - w - j >= 0 && !(&(c[i - w - j]) == nullptr)) {
                yansehuancun = yansehuancun + c[i - w - j];
                a++;
            }
            if (y < h - 1 - j && i + w + j <= w * h && !(&(c[i + w + j]) == nullptr)) {
                yansehuancun = yansehuancun + c[i + w + j];
                a++;
            }
            if (x > 0 + j && i - 1 - j <= w * h && !(&(c[i - 1 - j]) == nullptr)) {
                yansehuancun = yansehuancun + c[i - 1 - j];
                a++;
            }
            if (x < w - (1 + j) && i + 1 + j <= w * h && !(&(c[i + 1 + j]) == nullptr)) {
                yansehuancun = yansehuancun + c[i + 1 + j];
                //a++;
            }
        }
    }
    if (yansehuancun.x > 0 && yansehuancun.y > 0 && yansehuancun.z > 0)
        return c[i] = yansehuancun * (((1 - p) / 2.5 + p) / a * 1.5);
}
Vec DeNoisy2(int x, int y, int i, int w, int h, int samps, Vec c[], int i1=0) {
    Vec yansehuancun;
    double p = c[i].x > c[i].y && c[i].x > c[i].z ? c[i].x : c[i].y > c[i].z ? c[i].y : c[i].z;
    double a = 0;
    for (int j = 1; j <= i1; j++) {
        if (c[i].x <= jz && c[i].y <= jz && c[i].z <= jz && x > 0 && y > 0 && samps <= 64) {//降噪
            if (y > 0 + j + 1 && i - w - j - 1 <= w * h && i - w - j >= 0 && !(&(c[i - w - j]) == nullptr)) {
                yansehuancun = yansehuancun + c[i - w - j - 1];
                a++;
            }
            if (y < h - 1 - j - 1 && i + w + j + 1 <= w * h && !(&(c[i + w + j]) == nullptr)) {
                yansehuancun = yansehuancun + c[i + w + j + 1];
                a++;
            }
            if (x > 0 + j + w && i - 1 - j - w <= w * h && !(&(c[i - 1 - j]) == nullptr)) {
                yansehuancun = yansehuancun + c[i - 1 - j - w];
                a++;
            }
            if (x < w - (1 + j - 1) && i + 1 + j + w <= w * h && !(&(c[i + 1 + j]) == nullptr)) {
                yansehuancun = yansehuancun + c[i + 1 + j + w];
                //a++;
            }
        }
    }
    if (yansehuancun.x > 0 && yansehuancun.y > 0 && yansehuancun.z > 0)
        return c[i] = yansehuancun * (((1 - p) / 2.5 + p) / a * 1.5);
}
bool js( Vec a,Vec b) {
    if (fabs(a.x - b.x) > 0.02 && fabs(a.y - b.x) > 0.02 && fabs(a.z - b.z) > 0.02) {
        return true;
    }
    else {
        return false;
    }
}
void setRGBColor(Vec c) {
    int r_ = c.x * 255;
    int g_ = c.y * 255;
    int b_ = c.z * 255;
    cout << "\033[38;2;" << r_ << ";" << g_ << ";" << b_ << "m"; // 设置前景色为RGB(0-1)
}
//全屏显示
void MaxScreen() {
    HWND Hwnd = GetForegroundWindow();
    ShowWindow(Hwnd, SW_MAXIMIZE);
}
void SetPos(int x, int y)
{
    COORD point = { x, y };//光标要设置的位置x,y
    HANDLE HOutput = GetStdHandle(STD_OUTPUT_HANDLE);//使用GetStdHandle(STD_OUTPUT_HANDLE)来获取标准输出的句柄
    SetConsoleCursorPosition(HOutput, point);//设置光标位置
}
#define angle .1
int main(int argc, char* argv[]) {
    int r_;
    int g_;
    int b_;
    HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
    CONSOLE_SCREEN_BUFFER_INFO consoleInfo;

    // 获取当前控制台属性以便恢复
    GetConsoleScreenBufferInfo(hConsole, &consoleInfo);
    WORD originalAttrs = consoleInfo.wAttributes;

    // 设置控制台属性来实现RGB颜色
    SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_BLUE);
    MaxScreen();
    omp_set_num_threads(36);
    int w = 70, h = 40, samps = argc == 2 ? atoi(argv[1]) / 4 : 9, samps2 = samps; //图像大小及采样次数
    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // 摄像机位置和方向
    /*cx表示每次横向移动时移动多少,cy表示往下移动时移动多少,r为每次采样时记录的结果,c是用来计算最终图像的数组(↓)*/
    Vec cx = Vec(w * .5135 / h), cy = (cx % cam.d).norm() * .5135, r, * c = new Vec[w * h];

    std::ios::sync_with_stdio(false);
    // 创建字符缓冲区
    char* str = new char[w * h]; // 加上每行的换行符数量
    int ch;
    while (true) {
#pragma omp parallel for schedule(dynamic, 1) private(r)//多线程
        for (int y = 0; y < h; y++) { // 循环遍历图像的每一行
            //fprintf(stderr, "\r渲染中... %5.2f%%", 100. * y / (h - 1));//输出进度
            for (unsigned short x = 0, Xi[3] = { 0, 0, (unsigned short)(y * y * y) }; x < w; x++)
                for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++) //2x2子像素行
                    for (int sx = 0; sx < 2; sx++, r = Vec()) { //2x2子像素列,r用来记录本次获得颜色值
                        for (int s = 0; s < samps; s++) { //采样
                            /*滤波器,让采样位置在中间的概率大,在周围的概率小(↓)*/
                            double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                            /*r1,r2是0-2的数,有一半概率生成sqrt(r) - 1,一半生成1 - sqrt(2 - r)*/
                            double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                            //像素发出光线的方向(↓)
                            Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) + cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
                            /* 把摄像机的光线向前推以开始在内部, cam.o + d * 140是射线起点(相机和屏幕之间的距离是140个单位),
                            * (1. / samps)的意义是计算光照贡献的平均值,因为r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0, Xi)计算的是
                            所有采样结果的总和,否则光线亮度可能过大*/
                            r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0, Xi) * (1. / samps);
                            //x,y,索引,宽,高,采样次数,储存图像的数组,迭代次数
                            DeNoisy(x, y, i, w, h, samps, c, 1);
                            //x,y,索引,宽,高,采样次数,储存图像的数组,迭代次数
                            //DeNoisy2(x, y, i, w, h, samps, c, 1);
                        }
                        c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25; //因为是2*2的子像素,每个结果只占1/4,所以乘0.25
                        if (c[i - 1].x > 0.95 && c[i - 1].y > 0.95 && c[i - 1].z > 0.95) { //重要性采样
                            //samps = 1;
                            c[i - 1] = Vec(1, 1, 1);//Debug
                        }
                        else if (c[i - 1].x < 0.05 && c[i - 1].y < 0.05 && c[i - 1].z < 0.05) {
                            c[i - 1] = Vec();//Debug
                        }
                        else {
                            samps = samps2;
                        }
                    }
        }
            if (GetAsyncKeyState(VK_ESCAPE) & 0x8000) {
                return 0;
            } //当按下ESC时循环，ESC键的键值时27.
            else if (GetAsyncKeyState('A') & 0x8000) {
                cam.o.x = cam.o.x - 10;
            }
            else if (GetAsyncKeyState('D') & 0x8000) {
                cam.o.x = cam.o.x + 10;
            }
            else if (GetAsyncKeyState('W') & 0x8000) {
                cam.o.z = cam.o.z - 10;
            }
            else if (GetAsyncKeyState('S') & 0x8000) {
                cam.o.z = cam.o.z + 10;
            }
            else if (GetAsyncKeyState('Q') & 0x8000) {
                cam.o.y = cam.o.y + 10;
            }
            else if (GetAsyncKeyState('E') & 0x8000) {
                cam.o.y = cam.o.y - 10;
            }
            else if (GetAsyncKeyState(VK_UP) & 0x8000) {//上
                cam.d = VecRotate(cam.d,angle,Vec(1,0,0));
            }
            else if (GetAsyncKeyState(VK_DOWN) & 0x8000) {//下
                cam.d = VecRotate(cam.d, -angle, Vec(1, 0, 0));
            }
            else if (GetAsyncKeyState(VK_LEFT) & 0x8000) {//左
                cam.d = VecRotate(cam.d, angle, Vec(0, 1, 0));
            }
            else if (GetAsyncKeyState(VK_RIGHT) & 0x8000) {//右
                cam.d = VecRotate(cam.d, -angle, Vec(0, 1, 0));
            }
        std::string output;
        for (int i = 0; i < w * h; i++) {
            if ((i - 1) % w == 0) {
                output += "\n";
            }
            r_ = c[i].x * 255;
            g_ = c[i].y * 255;
            b_ = c[i].z * 255;
            string color = format("\033[38;2;{};{};{}m@ ", r_, g_, b_);
            output += color;
        }
        cout << output << flush;
        for (int clear = 0; clear < w * h; clear++) {
            c[clear] = Vec();
        }
        //system("cls");
        SetPos(0, 0);
    }
    /*FILE* f = fopen("image.ppm", "w");//文件写入
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w * h; i++)
        fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));*/
}//▓

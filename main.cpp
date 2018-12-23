
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
struct Vec {
    double x, y, z;
    Vec(){x=y=z=0;}
    Vec(double a, double b, double c){x=a, y=b, z=c;}
    Vec operator- (Vec v){return Vec(x-v.x,y-v.y,z-v.z);}
    Vec operator+ (Vec v){return Vec(x+v.x,y+v.y,z+v.z);}
    Vec operator* (double d){return Vec(x*d,y*d,z*d);}
    Vec operator/ (double d){return Vec(x/d,y/d,z/d);}
    Vec normalize(){double mg = sqrt(x*x+y*y+z*z); return Vec(x/mg, y/mg, z/mg); }
};

double dot(Vec v, Vec b){return (v.x*b.x + v.y*b.y + v.z*b.z); }

struct Ray{
    Vec o;  //Origin
    Vec d;  //direction
    Ray(Vec i, Vec j){o=i; d=j;}
};

struct Sphere{
    Vec c;      //center
    double r;   //radius
    
    Sphere(Vec i, double j){c=i,r=j;}
    Vec getNormal(Vec pi) {return Vec ((pi-c)/r);}
    bool intersect(Ray ray, double &t){
        Vec o = ray.o;
        Vec d = ray.d;
        Vec oc = o-c;
        double b = 2*dot(oc,d);
        double c = dot(oc,oc) - r*r;
        double disc = b*b-4*c;
        if(disc < 0) return false;
        else{
            disc = sqrt(disc);
            double t0 = -b-disc;
            double t1 = -b+disc;
        
            t =  (t0 < t1)? t0 : t1;
            return true;
            
        }
    }
};

struct Color{
    double r,g,b;
    Color(){r=g=b=0;}
    Color(double i, double j, double k){r=i,g=j,b=k;}
    Color operator* (double d){return Color(r*d,g*d,b*d);}
    Color operator+ (Color v){return Color((v.r+r)/2,(v.g+g)/2,(v.b+b)/2);}

    
};


float charToFloat(char* c)
{
	char parse = c[0];
	float number = 0;
	int afterVirgul = 0;
	int howAfterVirgul = 10;
	int i = 0;
	while(parse != '\0')
	{
	if(!afterVirgul && parse != '.')
	{
		switch (parse)
		{
		case '0' : number = number * 10; break;
		case '1' : number = number * 10 + 1; break;
		case '2' : number = number * 10 + 2; break;
		case '3' : number = number * 10 + 3; break;
		case '4' : number = number * 10 + 4; break;
		case '5' : number = number * 10 + 5; break;
		case '6' : number = number * 10 + 6; break;
		case '7' : number = number * 10 + 7; break;
		case '8' : number = number * 10 + 8; break;
		case '9' : number = number * 10 + 9; break;
		}		
	}
	else
		if(parse == '.')
		{
			afterVirgul = 1;
		}
		else
		{
		switch(parse)
			{	
			case '0' : howAfterVirgul *= 10; break;
			case '1' : number = number + 1./howAfterVirgul; howAfterVirgul *= 10; break;
			case '2' : number = number + 2./howAfterVirgul; howAfterVirgul *= 10; break;
			case '3' : number = number + 3./howAfterVirgul; howAfterVirgul *= 10; break;
			case '4' : number = number + 4./howAfterVirgul; howAfterVirgul *= 10; break;
			case '5' : number = number + 5./howAfterVirgul; howAfterVirgul *= 10; break;
			case '6' : number = number + 6./howAfterVirgul; howAfterVirgul *= 10; break;
			case '7' : number = number + 7./howAfterVirgul; howAfterVirgul *= 10; break;
			case '8' : number = number + 8./howAfterVirgul; howAfterVirgul *= 10; break;
			case '9' : number = number + 9./howAfterVirgul; howAfterVirgul *= 10; break;
			}
		}	
	i++;
	parse = c[i];
	}
	
	return number;

}

Color phong(Color pixelColor, float ka, float kd, float ks, Vec L, Vec N, Vec V)
{
	Color white(255, 255, 255);
    	Color black(0, 0, 0);
	Vec R = (N * (dot(N, L) * 2) ) - L;

	Color ambiant = white * ka;
	Color diffus  = white * kd * dot(L,N);
	Color Specular = white * ks * pow(dot(R, V), 22);
	pixelColor = black  + ambiant + diffus + Specular;

	return pixelColor; 
}

int main(int argc, char* argv[])
{
    const int W = 500;  // Image width
    const int H = 500;  // Image height
    
    ofstream out("phong.ppm");
    out << "P3\n" << W << '\n' << H << '\n' << "255\n";
    
    if(argc < 1)
	std::cout << "You should put parameters Ka, Kd and Ks" << std::endl;
    else     if(argc < 2)
    		std::cout << "You should complete parameters Kd and Ks" << std::endl; 
    	     else     if(argc < 3)
        		std::cout << "You should complete parameter Ks" << std::endl; 			    


    Color pixel_col[H][W];
    
    Color white (255, 255, 255);
    Color red (255, 0, 0);
    Color blue(0, 0, 255);
    Sphere sphere (Vec(W/2.0,H/2.0,200), 150);
    Sphere light (Vec(250, 250, 0), 1);
    //For each pixel
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++)
        {
            // send a ray through each pixel
            Ray ray(Vec(x, y, 0), Vec(0, 0, 1));
            double t = 20000;
            
            //check  for intersections
            if(sphere.intersect(ray, t)){
                
                //point of intersections
                Vec pi = ray.o + ray.d * t;
                
                //color the pixel
                Vec L = light.c - pi;
                Vec N = sphere.getNormal(pi);
                double dt = dot(L.normalize(), N.normalize());
                
		float ka = charToFloat(argv[1]);
		float kd = charToFloat(argv[2]);
		float ks = charToFloat(argv[3]);
		Vec Eye = pi - Vec(0, 0, 1);
		//std::cout <<			 		

   		pixel_col[y][x] = phong(blue, ka, kd, ks, L.normalize(), N.normalize(), Eye.normalize());
 
                //pixel_col[y][x] = red + (white * dt)* charToFloat(argv[1]);
                if(pixel_col[y][x].r < 0) pixel_col[y][x].r = 0;
                if(pixel_col[y][x].g < 0) pixel_col[y][x].g = 0;
                if(pixel_col[y][x].b < 0) pixel_col[y][x].b = 0;
                if(pixel_col[y][x].r > 255) pixel_col[y][x].r = 255;
                if(pixel_col[y][x].g > 255) pixel_col[y][x].g = 255;
                if(pixel_col[y][x].b > 255) pixel_col[y][x].b = 255;
		//std::cout <<pixel_col[y][x].r << endl;
            }
            out <<(int) pixel_col[y][x].r << endl;
            out <<(int) pixel_col[y][x].g << endl;
            out <<(int) pixel_col[y][x].b << endl;
        }
	//std::cout << argv[1] <<charToFloat(argv[1]) << std::endl;
    return 0;
}





#include<iostream>
#include<cstdlib>
#include<cstring>
#include<time.h>
#define VL 512
#define r 1
#define PL VL/8
#define DATA_TYPE double
#define n PL/(sizeof(DATA_TYPE))
#define CHECK
#define top ((DATA_TYPE)0.4)
#define mid ((DATA_TYPE)0.1)
#define down ((DATA_TYPE)0.2)
#define left ((DATA_TYPE)0.1)
#define right ((DATA_TYPE)0.2)
#define load(M,i,j,m) M[i*m+j]
//#define OUT
using namespace std;
//SME每行的bit数=VL=512bits=64Bytes,每行64Bytes=16个float=8个double,
//那么可以计算64/16=4组16*16的float外积,可以计算64/8=8组8*8的double外积,
//且物理上每组外积计算在SME上的排布是岔开的，有利于减少延迟等待，或者手动unroll循环展开

// FMOPA <ZAda>.S, <Pn>/M, <Pm>/M, <Zn>.S, <Zm>.S
void FMOPA_X(DATA_TYPE* ZA,DATA_TYPE* Zn,DATA_TYPE *Zm,int N){
    bool sub_op = false;
    int dim=n;

    for(int row=0;row<dim;row++){
        for(int col=0;col<dim;col++){
            ZA[row*(N+2*r)+col]+=Zn[row]*Zm[col];
        }    
    }
}
void FMOPA_Y(DATA_TYPE* ZA,DATA_TYPE* Zn,DATA_TYPE *Zm,int N){
    bool sub_op = false;
    int dim=n;

    for(int row=0;row<dim;row++){
        for(int col=0;col<dim;col++){
            ZA[row*(N+2*r)+col]+=Zn[row*(N+2*r)]*Zm[col];
        }    
    }
}
//2d5p stencil A[t][i][j] = 1*A[t-1][i][j-1] + 2*A[t-1][i][j+1] + 3*A[t-1][i][j] + 4*A[t-1][i-1][j] + 5*A[t-1][i+1][j];
int main(int argc,char* argv[])
{
    int N = atoi(argv[1]);
    int T = atoi(argv[2]);
	// int T = 1;
	// int N = 16;
	// int r = 1;
	//移动得到n+2r组长度为n的系数向量??(n+2r是不是可以减少操作)
    DATA_TYPE* h = (DATA_TYPE*)calloc((n+2*r)*n,sizeof(DATA_TYPE));
    DATA_TYPE* v = (DATA_TYPE*)calloc((n+2*r)*n,sizeof(DATA_TYPE));
    DATA_TYPE* A = (DATA_TYPE*)calloc((N+2*r)*(N+2*r),sizeof(DATA_TYPE));//源矩阵
    DATA_TYPE** dest = (DATA_TYPE**)malloc(sizeof(DATA_TYPE*)*2);//目标矩阵
    for(int i=0;i<2;i++)
    {
        *(dest+i) = (DATA_TYPE*)calloc((N+2*r)*(N+2*r),sizeof(DATA_TYPE));
    }
    int idx1=-2,idx2=-1,idx3=0;
    for(int i=0;i<n+2*r;i++)
    {
        if(i+idx1>=0 && i+idx1<n)load(h,i,i+idx1,n) = right;
        // if(i+idx2>=0 && i+idx2<n)h[i][i+idx2] = 3;
        if(i+idx3>=0 && i+idx3<n)load(h,i,i+idx3,n) = left;
    }
    for(int i=0;i<n+2*r;i++)
    {
        if(i+idx1>=0 && i+idx1<n)load(v,i,i+idx1,n) = down;
        if(i+idx2>=0 && i+idx2<n)load(v,i,i+idx2,n) = mid;
        if(i+idx3>=0 && i+idx3<n)load(v,i,i+idx3,n) = top;
    }
    srand((unsigned)time(NULL));
    for(int i=r;i<(N+r);i++)
    {    
        for(int j=r;j<(N+r);j++)
        {
            // A[i*(N+2*r)+j] = 1;
            load(A,i,j,(N+2*r)) = rand()%1024;
        }
    }
    cout << "矩阵A:" << endl;
    for(int i=0;i<N+2*r;i++)
	{
		for(int j=0;j<N+2*r;j++)
		{
			cout << load(A,i,j,(N+2*r)) << " ";
		}
		cout << endl;
	}
    cout << "----------------------------------" << endl;
#ifdef CHECK
    DATA_TYPE** comp = (DATA_TYPE**)malloc(sizeof(DATA_TYPE*)*2);
    for(int i=0;i<2;i++)
    {
        *(comp+i) = (DATA_TYPE*)calloc((N+2*r)*(N+2*r),sizeof(DATA_TYPE));
    }
    for(int i=r;i<N+r;i++)
    {
        for(int j=r;j<N+r;j++)
        {
            comp[0][i*(N+2*r)+j] =  top*A[(i-1)*(N+2*r)+j] + mid*A[i*(N+2*r)+j]+ down*A[(i+1)*(N+2*r)+j] + left*A[i*(N+2*r)+j-1] + right*A[i*(N+2*r)+j+1];
        }
    }
#ifdef OUT
    cout << "comp: " << endl;
    for(int i=0;i<N+2*r;i++)
	{
		for(int j=0;j<N+2*r;j++)
		{
			cout << comp[0][i*(N+2*r)+j] << " ";
		}
		cout << endl;
	}
    cout << "----------------------------------" << endl;
#endif
    for(int t=1;t<T;t++)
    {
        for(int i=r;i<N+r;i++)
        {
            for(int j=r;j<N+r;j++)
            {
                comp[t%2][i*(N+2*r)+j] =  top*comp[(t-1)%2][(i-1)*(N+2*r)+j] + mid*comp[(t-1)%2][i*(N+2*r)+j] + down*comp[(t-1)%2][(i+1)*(N+2*r)+j] + left*comp[(t-1)%2][i*(N+2*r)+j-1] + right*comp[(t-1)%2][i*(N+2*r)+j+1];
            }
        } 
#ifdef OUT
        cout << "comp: " << endl;
        for(int i=0;i<N+2*r;i++)
        {
            for(int j=0;j<N+2*r;j++)
            {
                cout << comp[t%2][i*(N+2*r)+j] << " ";
            }
            cout << endl;
        }
        cout << "----------------------------------" << endl;       
#endif
    }
#endif
	DATA_TYPE* temp = dest[0];
    for(int ii=0;ii<N;ii=ii+n)
    {
        dest[0] = temp + r*(N + 2*r) + r;; 
        dest[0] +=ii*(N+2*r);
        for(int jj=0;jj<N;jj=jj+n)
        {	
            for(int i = 0;i<n+2*r;i++)//是否此处要循环展开n次,带上边界需要n+2r组，sve不够用的重用问题
            {           
                //从A中取出向量MOV <Zd>.B, <Pg>/M, ZA0<HV>.B[<Ws>, <imm>]
                //计算系数向量v*格点横向量
                FMOPA_X(dest[0],v+i*n,A+(ii+i)*(N+2*r)+jj+r,N);
            }
            dest[0] += n; 				
        }
    }
    for(int ii=0;ii<N;ii=ii+n)
    {
        dest[0] = temp + r*(N + 2*r) + r;; 
        dest[0] +=ii*(N+2*r);
        for(int jj=0;jj<N;jj=jj+n)
        {	
            for(int i = 0;i<n+2*r;i++)//是否此处要循环展开n次,带上边界需要n+2r组，sve不够用的重用问题
            {           
                //从A中取出向量MOV <Zd>.B, <Pg>/M, ZA0<HV>.B[<Ws>, <imm>]
                //格点纵向量*计算系数向量h
                FMOPA_Y(dest[0],A+(ii+1)*(N+2*r)+jj+i,h+i*n,N);
            }
            dest[0] += n; 				
        }
    }
    dest[0] = temp;
#ifdef OUT
    cout << "dest: " << endl;
    for(int i=0;i<(N+2*r);i++)
	{
		for(int j=0;j<(N+2*r);j++)
		{
			cout << dest[0][i*(N+2*r)+j] << " ";
		}
		cout << endl;
	}
    cout << "---------------------------" << endl;
#endif
    for(int t=1;t<T;t++)
    {
        for(int i=0;i<(N+2*r);i++)
        {
            for(int j=0;j<(N+2*r);j++)
            {
                dest[(t)%2][i*(N+2*r)+j]  = 0;
            }
        }        
        temp = dest[t%2];
        for(int ii=0;ii<N;ii=ii+n)
        {
            dest[t%2] = temp + r*(N + 2*r) + r;; 
            dest[t%2] +=ii*(N+2*r);
            for(int jj=0;jj<N;jj=jj+n)
            {	
                for(int i = 0;i<n+2*r;i++)//是否此处要循环展开n次,带上边界需要n+2r组，sve不够用的重用问题
                {           
                    //从A中取出向量MOV <Zd>.B, <Pg>/M, ZA0<HV>.B[<Ws>, <imm>]
                    //计算系数向量v*格点横向量
                    FMOPA_X(dest[t%2],v+i*n,dest[(t-1)%2]+(ii+i)*(N+2*r)+jj+r,N);
                }
                dest[t%2] += n; 				
            }
        }
        for(int ii=0;ii<N;ii=ii+n)
        {
            dest[t%2] = temp + r*(N + 2*r) + r;; 
            dest[t%2] +=ii*(N+2*r);
            for(int jj=0;jj<N;jj=jj+n)
            {	
                for(int i = 0;i<n+2*r;i++)//是否此处要循环展开n次,带上边界需要n+2r组，sve不够用的重用问题
                {           
                    //从A中取出向量MOV <Zd>.B, <Pg>/M, ZA0<HV>.B[<Ws>, <imm>]
                    //格点纵向量*计算系数向量h
                    FMOPA_Y(dest[t%2],dest[(t-1)%2]+(ii+1)*(N+2*r)+jj+i,h+i*n,N);
                }
                dest[t%2] += n; 				
            }
        }

        dest[t%2] = temp;
#ifdef OUT
        cout << "dest: " << endl;
        for(int i=0;i<(N+2*r);i++)
        {
            for(int j=0;j<(N+2*r);j++)
            {
                cout << dest[t%2][i*(N+2*r)+j] << " ";
            }
            cout << endl;
        }
        cout << "---------------------------" << endl;
#endif
    }


#ifdef CHECK
    bool flag = 0;
    for(int i=0;i<(N+2*r);i++)
	{
		for(int j=0;j<(N+2*r);j++)
		{
            if(dest[(T-1)%2][i*(N+2*r)+j] != comp[(T-1)%2][i*(N+2*r)+j])
            {
                cout << "Wrong: " << dest[(T-1)%2][i*(N+2*r)+j] << "  Right: " << comp[(T-1)%2][i*(N+2*r)+j];
                flag = 1;
            }
		}
	}
    if(flag==0)cout << "Correct!" << endl;
#endif
    free(h);
    free(v);
    free(A);
    free(dest);
    free(comp);
    return 0;
}

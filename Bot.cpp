# Be made by Superbia_zyb

#pragma GCC optimize ("O2")
#pragma warning(disable:4996)
#include <bits/stdc++.h>
#include <immintrin.h>
#include "jsoncpp/json.h"
using namespace std;
#define _CMP_LT_OQ    0x11 
#define _CMP_GT_OQ    0x1e 
#define _CMP_EQ_OS    0x10

#define TIME (double)0.97*CLOCKS_PER_SEC*1

float f[300] __attribute__((aligned(32))) = {0};
double C = 0.6;
double CC =4000.6;
int total;
int dx[] = { 0,1,0,-1 }, dy[] = { -1,0,1,0 };
int mp[20][20],st_mp[20][20];
int l[2], r[2];
double D = 0.5;
int N, M;
int TT;
int ttp;
struct your
{
    int x,y;
    int col;
}snake[2][1000] ,q[1000];
struct my
{
    int x,y,dir;
}sp[1000][5];
int lp[400];
void look();
class game
{
public:
    int is(int rnd)  //本回合是否生长
    {
        if (rnd <= 9) return 1;
        if ((rnd - 9) % 3 == 0) return 1;
        return 0;
    }
    void bfs(your st, your ed,float vis[],float dis[])
    {
        int ll=1,rr=0;
        your tmp;
        tmp = st, tmp.col = -1,q[++rr]=tmp, vis[tmp.x * M + tmp.y] = 1, dis[tmp.x * M + tmp.y] = 0;
        tmp = ed, tmp.col = 1,q[++rr]=tmp, vis[tmp.x * M + tmp.y] = 1, dis[tmp.x * M + tmp.y] = 0;
        while (ll<=rr)
        {
            tmp = q[ll],ll++;
            int np = tmp.x * M + tmp.y;
            for (int dxx,dyy,dzz,i = 1; i <= lp[np] ; ++i)
            {
                dxx = sp[np][i].x, dyy = sp[np][i].y, dzz = dxx * M + dyy;
                if(mp[dxx][dyy]!=0) continue;
                if(vis[dzz]==0)
                {
                    vis[dzz] = tmp.col, dis[dzz] = dis[np] + 1.0;
                    your nmp; nmp.x = dxx, nmp.y = dyy, nmp.col = tmp.col, rr++, q[rr]=nmp;
                }
                else if(dis[dzz] == dis[np] + 1 && vis[dzz]!=tmp.col) vis[dzz]=2;
            }
        }
    }
    double S1(int col)
    {
        int aa[10];
        int myop = check(col,aa);
        int youop = check(-col,aa);

        int id = (col==-1)?1:-1;
        if(!myop && !youop) return -0.25*id;
        if(!myop) return -0.7;
        if(!youop) return 0.7;
        
        float vis[208] __attribute__((aligned(32))) = {0};
        float dis[208] __attribute__((aligned(32))) = {0};
        memset(dis, 0x4f,sizeof dis);
        bfs(snake[1][r[1]],snake[0][r[0]],vis,dis);
        
        __m256 al,ar,bmp,cmp,tmp,all_one,bp,cp,vll;

        bp = _mm256_set1_ps(-1);
        cp = _mm256_set1_ps(1);
        all_one = _mm256_set1_ps(1);
        al = _mm256_set1_ps(0);
        ar = _mm256_set1_ps(0);
        
        for(int i=0;i<TT;i+=8)
        {
            tmp = _mm256_load_ps(vis+i);
            vll = _mm256_load_ps(f+i);

            bmp = _mm256_cmp_ps(tmp,bp,_CMP_EQ_OS);
            bmp = _mm256_and_ps(bmp,all_one);
            bmp = _mm256_mul_ps(bmp,vll);
            
            cmp = _mm256_cmp_ps(tmp,cp,_CMP_EQ_OS);
            cmp = _mm256_and_ps(cmp,all_one);
            cmp = _mm256_mul_ps(cmp,vll);

            al = _mm256_add_ps(al,bmp);
            ar = _mm256_add_ps(ar,cmp);            
        }
        //look();
        float ans1=0.0,ans2=0.0;
        for(int i=0;i<8;i++)
        {
            ans1 = ans1 + al[i];
            ans2 = ans2 + ar[i];
        }
        return (ans1 - ans2)/(ans1 + ans2)*id;
    } 
    
    inline void del(int col)
    {
        int id = (col == -1);
        mp[snake[id][l[id]].x][snake[id][l[id]].y] = 0;
        l[id]++;
    }
    void move(int col, int dir)
    {
        int id = (col == -1);
        your &tmp = snake[id][r[id]];
        int dxx = dx[dir] + tmp.x, dyy = dy[dir] + tmp.y;
        r[id]++, snake[id][r[id]].x = dxx , snake[id][r[id]].y = dyy, mp[dxx][dyy] = col;
    }
    inline void remove(int col)
    {
        int id = (col == -1);
        your tmp = snake[id][r[id]];
        r[id]--, mp[tmp.x][tmp.y] = 0;
    }
    int check(int col, int a[])
    {
        int cnt = 0;
        int id = (col == -1);
        your &tmp = snake[id][r[id]];
        int np = tmp.x * M + tmp.y;
        for (register int i = 1; i <= lp[np] ; ++i)
        {
            if (mp[sp[np][i].x][sp[np][i].y] != 0) continue;
            a[++cnt] = sp[np][i].dir;
        }
        return cnt;
    }
}g;
void look()
{
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < M; ++j)
            printf("%2d ", mp[i][j]);
        printf("\n");
    }
    printf("\n");
}
class node
{
public:
    node() {}
    node(node* x, int tmp, int dxx, int dyy) { dir = tmp, col = -(x->col), fa = x, al = 0, rnd = x->rnd + 1; xx = dxx, yy = dyy; }
    node* son[90];
    node* fa = NULL;
    int n = 0, al = 0, dir, rnd, col;
    int xx, yy;
    double win = 0, vis = 0;
    double simulate(int depth, int rnd, int col)
    {
        int lp = -col;
        for (register int i = 1; i <= depth; ++i)
        {
            int op[10];
            int cnt = g.check(lp, op);
            if (cnt == 0)
            {
                if(rnd&1) g.remove(-lp);
                return g.S1(col);
            }
            int tmp = op[rand()%cnt + 1];
            g.move(lp, tmp);
            lp = -lp, rnd++;
            if (rnd % 2 == 0 && !g.is(rnd >> 1)) g.del(-1), g.del(1);
        }
        return g.S1(col);
    }
    void expand()
    {
        if (al) return;
        int id = (-col == -1);
        your tmp = snake[id][r[id]];
        int np = tmp.x * M + tmp.y;
        for (register int i = 1; i <= lp[np] ; ++i)
        {
            if (mp[sp[np][i].x][sp[np][i].y] != 0) continue;
            node* nmp = new node(this, sp[np][i].dir, sp[np][i].x,sp[np][i].y);
            son[++n] = nmp;
        }
        al = 1;
    }
    node* select()
    {
        node* tmp = this;
        if (tmp == NULL) return NULL;
        if (!tmp->vis)
        {
            if (tmp->rnd % 2 ==0)
            {
                if (!g.is(tmp->rnd >> 1))
                    g.del(-1), g.del(1);
            }
            return tmp;
        }
        while (tmp->n > 0)
        {
            tmp = tmp->get_best();
            g.move(tmp->col, tmp->dir);
            if (tmp->rnd % 2 == 0)
            {
                if (!g.is(tmp->rnd >> 1))
                    g.del(-1), g.del(1);
            }
        }
        return tmp;
    }
    node* get_best()
    {
        double tp = -100000.0, v;
        node* tmp = son[1];
        for (register int i = 1; i <= n; ++i)
        {
            if (son[i]->vis)
            {
                double nmp = 1.0/son[i]->vis;
                v = (son[i]->win * nmp) + C * sqrt(2 * log(vis) * nmp);
            }
            else v = 200000000;
            if (v > tp) tp = v, tmp = son[i];
        }
        return tmp;
    }
    node* get_ans()
    {
        if (!n) return NULL;
        node* tmp = son[1];
        for (int i = 1; i <= n; ++i)
            if (tmp->win * son[i]->vis < tmp->vis * son[i]->win) tmp = son[i];
        return tmp;
    }
    void update(double x)
    {
        node* tmp = this;
        while (tmp) tmp->vis++, tmp->win += x, x = -x, tmp = tmp->fa;
    }
}root;
void pi()
{
    for (int i = 1; i <= root.n; ++i)
    {
        if (!root.son[i]->vis) continue;
        printf("\n%d %d %lf %lf %lf %lf %lf", root.son[i]->xx, root.son[i]->yy, root.son[i]->win, root.son[i]->vis,
            (root.son[i]->win / root.son[i]->vis), C * sqrt(2 * log(root.vis) / root.son[i]->vis),
            (root.son[i]->win / root.son[i]->vis) + C * sqrt(2 * log(root.vis) / root.son[i]->vis));
    }
    printf("\n");
}
char de[1000000];
void work()
{
    if (!g.is(total)) g.del(-1), g.del(1);
    int len = (r[0] - l[0] + 1);
    for (register int i = l[0]; i <= r[0]; ++i)
        snake[0][i - l[0] + 1] = snake[0][i];
    for (register int i = l[1]; i <= r[1]; ++i)
        snake[1][i - l[1] + 1] = snake[1][i];
    r[0] = r[0] - l[0] + 1, l[0] = 1;
    r[1] = r[1] - l[1] + 1, l[1] = 1;

    memcpy(st_mp,mp,sizeof mp);
    int depth = 6;
    root.rnd = total << 1, root.col = 1, root.vis = 1, root.expand();
    if (root.n == 0)
    {
        printf("-1 -1");
        return;
    }
    double start = clock();
    for (int i=1;i<=root.n;i++)
    {
    	for(int j=1;j<=100;j++)
    	{
    	    g.move(root.son[i]->col, root.son[i]->dir);
    		node* leaf = root.son[i]->select();
	        if (leaf->vis)
	        {
	            leaf->expand();
	            if (leaf->n)
	                leaf = leaf->select();
	        }
	        double flag = leaf->simulate(depth + (leaf->rnd & 1), leaf->rnd, leaf->col);
	        leaf->update(flag);
	        memcpy(mp,st_mp,sizeof st_mp);
	        l[0] = l[1] = 1;
	        r[0] = r[1] = len;	
		}
	}
    for (; clock() - start < TIME;)
    {
        //pi();
        node* leaf = root.select();
        if (leaf->vis)
        {
            leaf->expand();
            if (leaf->n)
                leaf = leaf->select();
        }
        double flag = leaf->simulate(depth + (leaf->rnd & 1), leaf->rnd, leaf->col);
        leaf->update(flag);
        memcpy(mp,st_mp,sizeof st_mp);
        l[0] = l[1] = 1;
        r[0] = r[1] = len;
    }
    //pi();

    node* tmp = root.get_ans();
    Json::Value ret;
    ret["response"]["direction"] = tmp->dir;

    sprintf(de, "%lf ", root.vis);
    for (int i = 1; i <= root.n; ++i)
        sprintf(de + strlen(de), " ;%d %d %.0lf %.0lf",
            root.son[i]->xx, root.son[i]->yy, root.son[i]->win, root.son[i]->vis);

    ret["debug"] = de;
    Json::FastWriter writer;
    cout << writer.write(ret) << endl;
}
void input();
int main()
{
    // freopen("work.in", "r", stdin);
    // freopen("work.out", "w", stdout);
    input();
    work();
}
void input()
{

    srand((unsigned)time(0));

    l[0] = l[1] = 1;
    r[0] = r[1] = 0;

    string str, temp;
    getline(cin, temp);
    str += temp;
    Json::Reader reader; Json::Value input;
    reader.parse(str, input);
    M = input["requests"][(Json::Value::UInt) 0]["height"].asInt();  //棋盘高度
    N = input["requests"][(Json::Value::UInt) 0]["width"].asInt();   //棋盘宽度
    TT = ceil(N*M/8.0)*8;
    int x = input["requests"][(Json::Value::UInt) 0]["x"].asInt();  //读蛇初始化的信息
    if (x == 1)
    {
        r[1]++,snake[1][r[1]].x = snake[1][r[1]].y = 0, mp[0][0] = -1;
        r[0]++,snake[0][r[0]].x = N - 1, snake[0][r[0]].y = M-1, mp[N - 1][M - 1] = 1;
    }
    else
    {   
        r[0]++,snake[0][r[0]].x = snake[0][r[0]].y = 0, mp[0][0] = 1;
        r[1]++,snake[1][r[1]].x = N - 1, snake[1][r[1]].y = M-1, mp[N - 1][M - 1] = -1;
    }
    //处理地图中的障碍物
    int nmp = input["requests"][(Json::Value::UInt) 0]["obstacle"].size();
    for (int x, y, i = 0; i < nmp; ++i)
    {
        x = input["requests"][(Json::Value::UInt) 0]["obstacle"][(Json::Value::UInt) i]["x"].asInt();
        y = input["requests"][(Json::Value::UInt) 0]["obstacle"][(Json::Value::UInt) i]["y"].asInt();
        x--, y--, mp[y][x] = 2;
    }
    //根据历史信息恢复现场
    total = input["responses"].size();
    for (int dire, i = 0; i < total; ++i)
    {
        if (!g.is(i)) g.del(-1), g.del(1);
        dire = input["responses"][i]["direction"].asInt();
        g.move(-1, dire);
        dire = input["requests"][i + 1]["direction"].asInt();
        g.move(1, dire);
    }
	for(int ll=0;ll<=(N+1)/2 && ll<=(M+1)/2; ll++)
		for(int i=ll;i<N-ll;i++)
			for(int j=ll;j<M-ll;j++)
			{
				double tmp = (N+1)/2.0,nmp = (M+1)/2.0;
				f[i*M+j] = (min((N+1)/2,(M+1)/2)-(int)(tmp - ll))/CC+0.7;
				if(fabs(i+1-tmp)<=3 && fabs(j+1-nmp)<=3) f[i*M+j] += 0.1;
			}
    for(int i=0;i<N;i++)
        for(int cnt,j=0;j<M;j++)
        {
            cnt = 0;
            int tmp=i*M+j;
            for(int k=0;k<4;k++)
            {
                int x=dx[k]+i,y=dy[k]+j;
                if(x>=0 && x<N && y>=0 && y<M) 
                    cnt++,sp[tmp][cnt].x=x,sp[tmp][cnt].y=y,sp[tmp][cnt].dir=k;
            }
            lp[tmp]=cnt;
        }
}

<font size = 4>[Miller_Rabin](#1)</font>   
<font size = 4>[CRT & exgcd](#7)</font>  
<font size = 4>[Lucas](#8)</font>  
<font size = 4>[杜教筛](#14)</font>  
<font size = 4>[线性求逆元](#20)</font>

----------

<font size = 4>[树链剖分](#2)</font>   
<font size = 4>[LCA](#4)</font>  
&nbsp; &nbsp;  <font size = 3>[1.Tarjan](#5)</font>  
&nbsp; &nbsp;  <font size = 3>[2.倍增](#6)</font> 

<font size = 4>[Tarjan](#18)</font>  
&nbsp; &nbsp; <font size = 3>[1.强连通分量](#19)  
&nbsp; &nbsp; []()
</font>


<font size = 4>[Splay](#Splay)</font>

--------

<font size = 4>[可持久化线段树](#9)</font>  
&nbsp; &nbsp; <font size = 3>[静态第k小](#10)</font>  
&nbsp; &nbsp; <font size = 3>[动态第k小](#11)</font>

--------

<font size = 4>[KMP](#13)</font>  
<font size = 4>[最小表示法](#最小表示法)</font>  
<font size = 4>[Manacher](#Manacher)</font>  
<font size = 4>[ac-自动机](#12)</font>  
<font size = 4>[回文树](#16)</font>  
<font size = 4>[后缀数组](#17)</font>  
<font size = 4>[SAM](#SAM)</font>   
<font size = 4>[广义SAM](#GSAM)</font>
--------

<font size = 4>[FFT](#15)</font>

--------

<font size = 4>[快读、快写](#3)</font>  

----------


<h2 id = "1"> Miller_Rabin </h2>

```c++
#define ll unsigned long long
ll test[6] = {2, 3, 5, 233, 331};
ll qmul(ll a, ll b, ll mod){ // 快速模乘法
    ll ans = 0;
    while (b){
        if (b & 1){
            ans += a;
            ans %= mod;
        }
        a *= 2;
        a %= mod;
        b /= 2;
    }
    return ans;
}
ll qpow(ll a, ll n, ll mod)
{
    ll ret = 1;
    while (n)
    {
        if (n & 1)
            ret = qmul(ret, a, mod);
        a = qmul(a, a, mod);
        n >>= 1;
    }
    return ret;
}
bool Miller_Rabin(ll p)
{
    if (p < 2)
        return 0;
    if (p != 2 && p % 2 == 0)
        return 0;
    ll s = p - 1;
    while (!(s & 1))
        s >>= 1;
    for (int i = 0; i < 5; ++i)
    {
        if (p == test[i])
            return 1;
        ll t = s, m = qpow(test[i], s, p);
        while (t != p - 1 && m != 1 && m != p - 1)
        {
            m = qmul(m, m, p);
            t <<= 1;
        }
        if (m != p - 1 && !(t & 1))
            return 0;
    }
    return 1;
}
```
<h2 id = "7"> CRT & exgcd</h2>  
<font size = 3>  

$m_1, m_2....m_k$彼此互素，同余方程组$x=a_i(mod \ m_i)$在模$M=m_1m_2...m_k$下有为一解$x=a_1M_1M_1^{-1}+a_2M_2M_2^{-1}+...+a_kM_kM_k^{-1} \ (mod \ M)$  
$M_i=M/m_i, M_i^{-1}$是$M_i$模$m_i$的逆元

</font>

```c++
void extend_Euclid(int a, int b, int &x, int &y)
{
    if(b == 0){
        x = 1;
        y = 0;
        return;
    }
    extend_Euclid(b, a % b, x, y);
    int tmp = x;
    x = y;
    y = tmp - (a / b) * y;
}
int CRT(int a[],int m[],int n)
{
    int M = 1;
    int ans = 0;
    for(int i=1; i<=n; i++)
        M *= m[i];
    for(int i=1; i<=n; i++)
    {
        int x, y;
        int Mi = M / m[i];
        extend_Euclid(Mi, m[i], x, y);
        ans = (ans + Mi * x * a[i]) % M;
    }
    if(ans < 0) ans += M;
    return ans;
}
```

<h2 id = "8"> Lucus </h2>  

```c++
LL lucas(int n, int m) 
{
    if(n<m) return 0;
    LL ans = 1;
    for(; m; n/=P, m/=P) ans = ans *C (n%P, m%P) % P;
    // C(n, m) 组合数
    return ans;
}
```

<h2 id = 14> 杜教筛 </h2>

$ans_1=\sum_{i=1}^n\phi(i)$  
$ans_2=\sum_{i=1}^n\mu(i)$

```c++
map<LL, LL> _phi, _mu; // 考虑使用数组优化
LL phi[MAX+10], mu[MAX+10];
bool isprime[MAX+10];
vector<LL> prime;
void init()
{
    mset(isprime, 1);
    isprime[0] = isprime[1] = 0;
    phi[1] = 1;
    mu[1] = 1;
    for(int i = 2; i <= MAX; i++){
        if(isprime[i]){
            phi[i] = i-1;
            mu[i] = -1;
            prime.push_back(i);
        }
        for(int j = 0; j < prime.size() && i * prime[j] <= MAX; j++){
            isprime[i * prime[j]] = 0;
            if(i % prime[j] == 0){
                mu[i * prime[j]] = 0;
                phi[i * prime[j]] = phi[i] * prime[j];
                break;
            }
            mu[i * prime[j]] = -mu[i];
            phi[i * prime[j]] = phi[i] * (prime[j] - 1);
        }
    }
    for(int i = 2; i <= MAX; i++){
        mu[i] += mu[i-1];
        phi[i] += phi[i-1];
    }
}
LL S1(LL n)
{
    // phi
    if(n <= MAX){
        return phi[n];
    }
    map<LL, LL>::iterator it = _phi.find(n);
    if(it != _phi.end()){
        return it->second;
    }
    LL ans = (n + 1) * n >> 1;
    for(int i = 2, j; i <= n; i = j + 1){
        j = n / (n / i);
        ans -= (j - i + 1) * S1(n / i);
    }
    return _phi[n] = ans;
}
LL S2(LL n)
{
    // mu
    if(n <= MAX){
        return mu[n];
    }
    map<LL, LL>::iterator it = _mu.find(n);
    if(it != _mu.end()){
        return it->second;
    }
    LL ans = 1;
    for(int i = 2, j; i <= n; i = j+ 1){
        j = n / (n / i);
        ans -= (j - i + 1) * S2(n / i);
    }
    return _mu[n] = ans;
}
int main()
{
    init();
    LL T, n;
    scanf("%lld", &T);
    while(T--){
        scanf("%lld", &n);
        printf("%lld %lld\n", S1(n), S2(n));
    }
    return 0;
}
```

<h2 id = 20> 线性求逆元 </h2>

$$
令p=ki+r, k=\lfloor \frac{p}{i} \rfloor, r=p \ \% \ i \\
\begin{aligned}
    ki+r &=0\quad(mod \ \ p) \\
    ki &= -r \quad(mod \ \ p)\\
    i^{-1} &= (p - \lfloor \frac{p}{i} \rfloor)  * (p\%i)^{-1} \quad (mod\ \ p)
\end{aligned}

$$  

---------

<h2 id = 2> 树链剖分 </h2>

```C++  
using namespace std;
const int inf = 0x3f3f3f3f;
const LL __64inf = 0x3f3f3f3f3f3f3f3f;
const int MAX = 1e5 + 505;
const int Mod = 9999991;

int p;
struct edge
{
    int v, next;
}E[MAX<<1];
int tot, head[MAX];
void addedge(int u, int v){
    E[tot].v = v;
    E[tot].next = head[u];
    head[u] = tot++;
}
int w[MAX], wt[MAX];
struct node
{
    int val, lazy;
}st[MAX<<2];
void pushup(int rt){
    st[rt].val = (st[rt<<1].val + st[rt<<1|1].val) % p;
}
void pushdown(int rt, int l, int r){
    if(st[rt].lazy){
        int mid = (l + r) >> 1;
        st[rt<<1].val = (st[rt<<1].val + (mid-l+1) * st[rt].lazy) % p;
        st[rt<<1|1].val = (st[rt<<1|1].val + (r - mid) * st[rt].lazy) % p;
        // st[rt<<1|1].lazy = st[rt<<1].lazy = st[rt].lazy;
        (st[rt<<1|1].lazy += st[rt].lazy) % p;
        (st[rt<<1].lazy += st[rt].lazy) % p;
        st[rt].lazy = 0;
    }
}
void build(int rt, int l, int r){
    if(l == r){
        st[rt].val = wt[l] % p;
        return ;
    }
    int mid = (l + r) >> 1;
    build(rt<<1, l, mid);
    build(rt<<1|1, mid+1, r);
    pushup(rt);
}
void updata(int rt, int l, int r, int L, int R, int x){
    if(l >= L and r <= R){
        st[rt].lazy += x;
        st[rt].val = (st[rt].val + x * (r - l + 1)) % p;
        return ;
    }
    int mid = (l + r) >> 1;
    pushdown(rt, l, r);
    if(L <= mid)    updata(rt<<1, l, mid, L, R, x);
    if(R > mid)     updata(rt<<1|1, mid+1, r, L, R, x);
    pushup(rt);
}
int query(int rt, int l, int r, int L, int R){
    if(l >= L and r <= R){
        return st[rt].val ;
    }
    int mid = (l + r) >> 1, res = 0;
    pushdown(rt, l, r);
    if(L <= mid)    res += query(rt<<1, l, mid, L, R);
    if(R > mid)     res += query(rt<<1|1, mid+1, r, L, R);
    return res % p;
}

int tim;
int in[MAX], f[MAX], top[MAX], siz[MAX], dep[MAX], son[MAX];
void dfs1(int u, int p, int d){
    dep[u] = d;
    siz[u] = 1;
    f[u] = p;
    int Max = 0;
    for(int i = head[u]; ~i; i = E[i].next){
        int v = E[i].v;
        if(v == p)  continue;
        dfs1(v, u, d+1);
        if(siz[v] > Max){
            Max = siz[v];
            son[u] = v;
        }
        siz[u] += siz[v];
    }
}
void dfs2(int u, int p, int topf){
    in[u] = ++tim;
    top[u] = topf;
    wt[tim] = w[u];
    if(!son[u]) return ;
    dfs2(son[u], u, topf);
    for(int i = head[u]; ~i ; i = E[i].next){
        int v = E[i].v;
        if(v == son[u] || v == p)   continue;
        dfs2(v, u, v);
    }
}
int main()
{
    int n, m, r;
    mset(head, -1);
    // freopen("in.txt", "r", stdin);
    scanf("%d%d%d%d", &n, &m, &r, &p);
    rep(i, 1, n+1)
        scanf("%d", &w[i]);
    rep(i, 1, n){
        int u, v;
        scanf("%d%d", &u, &v);
        addedge(u, v);
        addedge(v, u);
    }
    dfs1(r, -1, 1);
    dfs2(r, -1, r);
    build(1, 1, n);
    rep(i, 0, m){
        int cas, x, y, z;
        scanf("%d", &cas);
        if(cas == 1){
            scanf("%d%d%d", &x, &y, &z);
            z %= p;
            while(top[x] != top[y]){
                if(dep[top[x]] > dep[top[y]]){
                updata(1, 1, n, in[top[x]], in[x], z);
                //    debug(st[r].val);
                x = f[top[x]];

                }
                else{
                    updata(1, 1, n, in[top[y]], in[y],z);
                    // debug(st[r].val);
                    y = f[top[y]];
                }
            }
            if(dep[x] > dep[y]) swap(x, y);
            updata(1, 1, n, in[x], in[y], z);
            // debug(st[r].val);
        }
        else if(cas == 2){
            scanf("%d%d", &x, &y);
            int res = 0;
            while(top[x] != top[y]){
                if(dep[top[x]] > dep[top[y]]){
                    res = (res + query(1, 1, n, in[top[x]], in[x])) % p;
                    x = f[top[x]];
                }
                else{
                    res = (res + query(1, 1, n, in[top[y]], in[y])) % p;
                    y = f[top[y]];
                }
            }
            if(dep[x] > dep[y]) swap(x, y);
            res = (res + query(1, 1, n, in[x], in[y])) % p;
            if(res < 0) res += p;
            printf("%d\n", res);
        }
        else if(cas == 3){
            scanf("%d%d", &x, &z);
            z %= p;
            // debug(query(1, 1, n, in[x], in[x] + siz[x] - 1));
            updata(1, 1, n, in[x], in[x] + siz[x] - 1, z);
            // debug(query(1, 1, n, in[x], in[x]+siz[x]-1));
        }   
        else{
            scanf("%d", &x);
            int res = query(1, 1, n, in[x], in[x]+siz[x]-1) % p;
            printf("%d\n", res);
        }
    }
    return 0;
}  
```

<h2 id = "4"> LCA </h2>  
<h3 id = "5"> 1.Tarjan </h3>
   
```c++
struct edge{
    int v, val, next;
}E[MAX<<1];
int head[MAX];
int tot;
void addedge(int u, int v, int val){
    E[tot].v = v;
    E[tot].val = val;
    E[tot].next = head[u];
    head[u] = tot++;
}
int q[MAX], ans[505];
void addquery(int u, int v, int id){
    E[tot].v = v;
    E[tot].val = id;
    E[tot].next = q[u];
    q[u] = tot++;
}
int f[MAX], dep[MAX];
bool vis[MAX];
int Find(int x){
    return x == f[x] ? x : f[x] = Find(f[x]);
}
void dfs(int u, int p, int d){
    dep[u] = d;
    for(int i = head[u]; ~i ; i = E[i].next){
        int v = E[i].v;
        if(v == p)  continue;
        dfs(v, u, d+E[i].val);
        f[Find(v)] = u;
    }
    vis[u] = 1;
    for(int i = q[u]; ~i; i = E[i].next){
        int v = E[i].v;
        if(vis[v] and ans[E[i].val] == -1)
            ans[E[i].val] = dep[v] + dep[u] - 2 * dep[Find(v)];
    }
}
int main()
{
    int n, m;
    int cas;
    scanf("%d", &cas);
    while (cas--)
    {
        scanf("%d%d", &n, &m);
        mset(ans, -1);
        mset(vis, 0);
        mset(head, -1);
        mset(q, -1);
        rep(i, 0, n+1) f[i] = i;
        tot = 0;
        rep(i, 1, n){
            int u, v, val;
            scanf("%d%d%d", &u, &v, &val);
            addedge(u, v, val);
            addedge(v, u, val);
        }
        rep(i, 0, m){
            int u, v;
            scanf("%d%d", &u, &v);
            addquery(u, v, i+1);
            addquery(v, u, i+1);
        }
        dfs(1, -1, 0);
        for(int i = 1; i <= m; i++)
            printf("%d\n", ans[i]);
    }
   return 0; 
}
```

<h3 id = "6"> 2.倍增 </h3>
   
```c++
struct edge
{
    int v, val, next;
} E[MAX << 1];
int head[MAX];
int tot;
void addedge(int u, int v, int val)
{
    E[tot].v = v;
    E[tot].val = val;
    E[tot].next = head[u];
    head[u] = tot++;
}
int dis[MAX], dep[MAX];
int f[MAX][20];
void dfs(int u, int p, int val, int d)
{
    f[u][0] = p;
    for(int i = 1; i < 20; i++)
        f[u][i] = f[f[u][i-1]][i-1];
    dis[u] = val, dep[u] = d;
    for(int i = head[u]; ~i; i = E[i].next){
        int v = E[i].v;
        if(v == p)  continue;
        dfs(v, u, val+E[i].val, d+1);
    }
}
int LCA(int u, int v){
    int res = dis[u] + dis[v];
    if(dep[u] < dep[v])
        swap(u, v);
    for(int i = 19; i >= 0; i--){
        if(dep[f[u][i]] >= dep[v])
            u = f[u][i];
    }
    if(u == v){
        res -= 2 * dis[v];
        return res;
    }
    for(int i = 19; i >= 0; i--){
        if(f[u][i] != f[v][i])
        {
            u = f[u][i];
            v = f[v][i];
        }
    }
    return res - 2 * dis[f[u][0]];
}
int main()
{
    int n, m;
    int cas;
    scanf("%d", &cas);
    while (cas--)
    {
        scanf("%d%d", &n, &m);
        mset(head, -1);
        tot = 0;
        rep(i, 1, n)
        {
            int u, v, val;
            scanf("%d%d%d", &u, &v, &val);
            addedge(u, v, val);
            addedge(v, u, val);
        }
        dfs(1, 0, 0, 1);
        rep(i, 0, m)
        {
            int u, v;
            scanf("%d%d", &u, &v);
            int ans = LCA(u, v);
            printf("%d\n", ans);
        }
    }
    return 0;
}
```

1.  ### RMQ ###  
   

<h2 id = 18> Tarjan </h2>
<h3 id = 19> 强连通分量</h3>

```c++
int n, m, tim;
vector<int> E[MAX];
stack<int> stk;
bool vis[MAX];
int dfn[MAX], low[MAX];
int col[MAX], num;
void tarjan(int u){
    tim++;
    dfn[u] = low[u] = tim;
    stk.push(u);
    vis[u] = 1;
    for(int i = 0; i < E[u].size(); i++){
        int v = E[u][i];
        if(!dfn[v]){
             tarjan(v);
            low[u] = min(low[u], low[v]);
        }
        else {
            if(vis[v])
                low[u] = min(low[u], low[v]);
        }
    }
    if(low[u] == dfn[u]){
        col[u] = ++num;
        vis[u] = 0;
        while(!stk.empty() and stk.top() != u){
            col[stk.top()] = num;
            vis[stk.top()] = 0;
            stk.pop();
        }
        stk.pop();
    }
}
int main()
{
    scanf("%d%d", &n, &m);
    for(int i = 1; i <= m; i++){
        int u, v;
        scanf("%d%d", &u, &v);
        E[u].push_back(v);
    }
    for(int i = 1; i <= n; i++){
        if(!dfn[i]) tarjan(i);
    }
    map<int, int> cnt;
    for(int i = 1; i <= n; i++){
        cnt[col[i]]++;
    }
    int ans = 0;
    for(auto it = cnt.begin(); it != cnt.end(); it++)   ans += it->second > 1;
    printf("%d\n", ans);
    return 0;
}
```

## Splay 
```c++
template<class T>
struct node{
    T val;
    int  l, r, p;
    node(int val = 0, int l = -1, int r = -1, int p=-1) : val(val), l(l), r(r), p(p){}
};
template<class T>
struct SplayTree{
    node<T> tr[MAX];
    int  rt, tot;;
    void init(){
        tot = 0, rt = -1;
    }
    void Splay(int u);
    void right_rotate(int u);
    void left_rotate(int u);
    int find(T x);
    void insert(T x);
    void erase(T x);
    int successor(int u);
    int predecessor(int u);
    node<T>& operator [](int x){
        return tr[x];
    }
};
template<class T>
void SplayTree<T>::Splay(int u){
    if(u == rt) return;
    for(int& p = tr[u].p; tr[u].p != -1; ){
        if(u == tr[p].l)
            right_rotate(p);
        else
            left_rotate(p);
    }
    rt = u;
}
template<class T>
void SplayTree<T>::right_rotate(int u){
    // u must have left child
    int v = tr[u].l;
    tr[u].l = tr[v].r;
    if(tr[v].r != -1) tr[tr[v].r].p = u;
    tr[v].r = u;
    if(tr[u].p != -1){
        int p = tr[u].p;
        if(tr[p].l == u) tr[p].l = v;
        else tr[p].r = v;
    }
    tr[v].p = tr[u].p;
    tr[u].p = v;
}
template<class T>
void SplayTree<T>::left_rotate(int u){
    // u must have right child
    int v = tr[u].r;
    tr[u].r = tr[v].l;
    if(tr[v].l != -1) tr[tr[v].l].p = u;
    tr[v].l = u;
    if(tr[u].p != -1){
        int p = tr[u].p;
        if(tr[p].l == u) tr[p].l = v;
        else tr[p].r = v;
    }
    tr[v].p = tr[u].p;
    tr[u].p = v;
}
template<class T>
int SplayTree<T>::find(T x){
    int u = rt;
    while(~u and tr[u].val != x){
        if(tr[u].val > x)
            u = tr[u].l;
        else u = tr[u].r;
    }
    if(u != -1) Splay(u);
    return u;
}
template<class T>
void SplayTree<T>::insert(T x){
    if(rt == -1){
        rt = 0;
        tr[tot++] = node<T>(x, -1, -1, -1);
        return;
    }
    int u = rt, p = -1;
    while(~u and tr[u].val != x){
        p = u;
        if(tr[u].val > x)
           u = tr[u].l;
        else u = tr[u].r;
    }
    if(u != -1) return ;
    tr[tot++] = node<T>(x, -1, -1, p);
    if(x < tr[p].val) tr[p].l = tot-1;
    else tr[p].r = tot-1;
    tr[tot-1].p = p;
    Splay(tot-1);
}
template<class T>
int SplayTree<T>::successor(int u){
    if(tr[u].r != -1){
        int v = tr[u].r;
        while(tr[v].l != -1)
            v = tr[v].l;
        return v;
    }
    else{
        int p = tr[u].p;
        while(p != -1 and u == tr[p].r) u = p, p = tr[u].p;
        return p;
    }
}
template<class T>
int SplayTree<T>::predecessor(int u){
    if(tr[u].l != -1){
        int v = tr[u].l;
        while(tr[v].r != -1) v = tr[v].r;
        return v;
    }
    else{
        int p =tr[u].p;
        while(p != -1 and u == tr[p].l) 
            u = p, p = tr[u].p;
        return p;
    }
}
template<class T>
void SplayTree<T>::erase(T x){
    int u = find(x);
    if(tr[rt].l != -1){
        tr[tr[rt].l].p = -1;
        int v = tr[u].l;
        while(tr[v].r != -1) v = tr[v].r;
        Splay(v);
        rt = v;
        tr[v].r = tr[u].r;
        if(tr[u].r != -1) tr[tr[u].r].p = v;
    }
    else if(tr[rt].r != -1) {
        tr[tr[rt].r].p = -1;
        int v = tr[u].r;
        while(tr[v].l != -1) v = tr[v].l;
        Splay(v);
        rt = v;
        tr[v].l = tr[u].l;
        if(tr[u].l != -1) tr[tr[u].l].p = v;
    } 
    else init();
}
```


---------

<h2 id = 9>可持久化线段树</h2>  
<h3 id = 10> 静态第k小 </h3>  

```c++  
struct node
{
    int ls, rs;
	int val;
};
node tr[MAX * 20];
int Rt[MAX];
int cur;
void pushup(int rt){
    // rt->val = rt->rs->val + rt->ls->val;
	tr[rt].val = tr[tr[rt].ls].val + tr[tr[rt].rs].val;
}
int build(int l, int r)
{
    int k = cur++;
    if(l == r){
        tr[k].val = 0;
        // rt->ls = rt->rs = NULL;
        return k;
    }
    int mid = (l + r) >> 1;
    tr[k].ls= build(l, mid);
    tr[k].rs = build(mid+1, r);
    pushup(k);
    return k;
}
int updata(int o, int l, int r, int idx, int val)
{
    int k = cur++;
	tr[k] = tr[o];
    if(l == r){
       	tr[k].val += val;
        return k;
    }
    int mid = (l + r) >> 1;
    if(idx <= mid)
        tr[k].ls = updata(tr[o].ls, l, mid, idx, val);
    else tr[k].rs = updata(tr[o].rs, mid+1, r, idx, val);
    pushup(k);
    return k;
}
int query(int rt1, int rt2, int l, int r, int kth)
{
    if(l == r){
        return l;
    }
    int mid = (l + r) >> 1;
    // int tmp = rt2->ls->val - rt1->ls->val ;
	int tmp = tr[tr[rt2].ls].val - tr[tr[rt1].ls].val;
    if(tmp >= kth){
        return query(tr[rt1].ls, tr[rt2].ls, l, mid, kth);
    }
    else
        return query(tr[rt1].rs, tr[rt2].rs, mid+1, r, kth-tmp);
}
int a[MAX], A[MAX];
int main()
{
    int n, m;
    while(~scanf("%d %d", &n, &m)){
		cur = 0;
    for(int i = 1; i <= n; i++){
        scanf("%d", &A[i]);
        a[i] = A[i];
    }
    sort(A+1, A+n+1);
    int cnt = 1;
    for(int i = 2; i <= n; i++){ 
        if(A[cnt] != A[i])
            A[++cnt] = A[i];
    }
    Rt[0] = build(1, cnt);
    for(int i = 1; i <= n; i++){
        int p = lower_bound(A+1, A+cnt+1, a[i]) - A;
        Rt[i] = updata(Rt[i-1], 1, cnt, p, 1);
    }
    for(int i = 1; i <= m; i++){
        int x, y, k;
        scanf("%d %d %d", &x, &y, &k);
        int id = query(Rt[x-1], Rt[y], 1, cnt, k);
        printf("%d\n", A[id]);
    }
	}
    return 0;
}
```

<h3 id = 11>动态第k小</h2>

```c++
int n, m, N;
int a[MAX];
int num[MAX<<1];
struct node
{
    int l, r,  val;
}st[MAX * 400];
struct op
{
    char c;
    int x, y, z;
};
vector<op> ops;
int bit[MAX], tot;
// vector<int> cnt[2];
int cnt[2], tmp[2][50];
void updata(int &rt, int l, int r, int pos, int x){
    if(!rt) rt = ++tot;
    st[rt].val += x;
    if(l == r){
        return;
    }
    int mid = (l + r) >> 1;
    if(pos <= mid) updata(st[rt].l, l, mid, pos, x);
    else updata(st[rt].r, mid+1, r, pos, x);
    // pushup(rt);
}
int query(int l, int r, int k){
    if(l == r){
        return l;
    }
    int mid = (l + r) >> 1, sum = 0;
    for(int i = 1; i <= cnt[1]; i++)  
        sum += st[st[tmp[1][i]].l].val;
    for (int i = 1; i <= cnt[0]; i++)
        sum -= st[st[tmp[0][i]].l].val;
    if(k <= sum) {
        for (int i = 1; i <= cnt[1]; i++)
            tmp[1][i] = st[tmp[1][i]].l;
        for (int i = 1; i <= cnt[0]; i++)
            tmp[0][i] = st[tmp[0][i]].l;
        return query(l, mid, k);
    }
    else{
        for(int i = 1; i <= cnt[1]; i++)
            tmp[1][i] = st[tmp[1][i]].r;
        for(int i = 1; i <= cnt[0]; i++)
            tmp[0][i] = st[tmp[0][i]].r;
        return query(mid+1, r, k-sum);
    }
}   
void add(int pos, int x){
    int k = lower_bound(num+1, num+1+N, a[pos])-num;
    for(int i = pos; i <= n; i += i & -i){
        updata(bit[i], 1, N, k, x);
    }
}
int Query(int l, int r, int k){
    // cnt[0].clear(); cnt[1].clear();
    mset(tmp, 0);
    cnt[0] = cnt[1] = 0;
    for(int i = r; i; i -= i & -i)
        tmp[1][++cnt[1]] = bit[i];
    for(int i = l-1; i; i -= i & -i)
        tmp[0][++cnt[0]] = bit[i];
    return query(1, N, k);
}
int main()
{
#ifdef DEBUG
    freopen("in.txt", "r", stdin);
#endif  
    scanf("%d%d", &n, &m);
    for(int i = 1; i <= n; i++){
        scanf("%d", a+i);
        num[++N] = a[i];
    }
    for(int i = 1; i <= m; i++){
        // char c = getchar();
        // while(c != 'C' and c != 'Q')    c = getchar();
        char c;
        cin >> c;
        if(c == 'C'){
            int x, y;
            scanf("%d%d", &x, &y);
            ops.push_back({c, x, y, 0});
            num[++N] = y;
        }
        else{
            int x, y, z;
            scanf("%d%d%d", &x, &y, &z);
            ops.push_back({c, x, y, z,});
        }
    }
    sort(num+1, num+1+N);
    N = unique(num+1, num+1+N)-num-1;
    for(int i = 1; i <= n; i++)
        add(i,  1);
    for(int i = 0; i < ops.size(); i++){
        if(ops[i].c == 'Q'){
            printf("%d\n", num[Query(ops[i].x, ops[i].y, ops[i].z)]);
        }
        else{
            add(ops[i].x, -1);
            a[ops[i].x] = ops[i].y;
            add(ops[i].x, 1);
        }
    }
    return 0;
}
```

----------

<h2 id = 13> KMP </h2>

```c++
string s, t;
int ans;
int nex[MAX];
void calc(){
    nex[0] = -1;
    int k = -1;
    for(int i = 1; i < t.size(); i++){
        while(k > -1 and t[k+1] != t[i])    
            k = nex[k];
        if(t[k+1] == t[i]) k++;
        nex[i] = k;
    }
}
void solve(){
    int k = -1;
    calc();
    for(int i = 0; i < s.size(); i++){
        while(k > -1 and t[k+1] != s[i])
            k = nex[k];
        if(t[k+1] == s[i])  k++;
        if(k == t.size()-1){
            ans++;
            k = nex[k];
        } 
    }
}
int main()
{
    int T;
    //freopen("in.txt", "r", stdin);
    ios::sync_with_stdio(false);
    cin >> T;
    while(T--){
        cin >> t >> s;
        ans = 0;
        solve();
        cout << ans << "\n";
    }
    return 0;
}
```

## 最小表示法 

```cpp

struct Less{
    bool operator () (char a, char b) {
        return a < b;
    }
};
struct Greator{
    bool operator () (char a, char b) {
        return a > b;
    }
};
template<typename Cmp>
int calcu(const string &s){
    int k = 0, i = 0, j = 1;
    int n = s.size();
    Cmp cmp;
    while (k < n and i < n and j < n)
    {
        if(s[(i+k)%n] == s[(j+k)%n]) k++;
        else {
            // s[(i+k)%n] < s[(j+k)%n] ? j = j + k + 1 : i = i + k + 1;
            cmp(s[(i+k)%n], s[(j+k)%n]) ? j = j + k + 1 : i = i + k + 1;
            if(i == j) i++;
            k = 0;
        }
    }
    return min(i, j);
}
```

## Manacher

> version 1
```c++
char s[MAX];
char t[(MAX<<1)+50];
int p[(MAX<<1)+50];
void init(){
	int n = strlen(s);
	t[0] = '@';
	t[1] = '#';
	t[2*n+2] = 0;
	for(int i = 0; s[i]; i++){
		t[(i<<1)+3] = '#';
		t[(i<<1)+2] = s[i];
	}
}
void solve(){
	init();
	// mset(p, 0);
	int Max = 0, id, ans = 0, ans_id;
	for(int i = 1; t[i]; i++){
		if(Max > i){
			int j = (id << 1) - i;
			p[i] = min(p[j], Max - i);
		}
		else p[i] = 1;
		for(; t[i-p[i]] == t[i+p[i]]; p[i]++);
		if(p[i]+i > Max) Max = p[i] + i, id = i;
		// ans = max(ans, p[i]);
		if(p[i] > ans) ans = p[i], ans_id = i;
	}
	printf("%d\n", --ans);
}
int main(){
#ifdef DEBUG
	freopen64("in", "r", stdin);
#endif
	int T = 1;
	while(~scanf("%s", s)){
		if(!strcmp(s, "END")) break;
		printf("Case %d: ", T++);
		solve();
	}
	return 0;
}
```

> version 2
```c++ 
int radius[MAX];
string manacherStr(const string &s){
    string ans = "#";
    for(int i = 0; i < s.size(); i++){
        ans += s[i];
        ans += "#";
    }
    return ans; 
}
int manacher(const string &s){
    if(s.size() == 0) return 0;
    int R = -1;
    int C = -1;
    int Max = -inf;
    string str = manacherStr(s);
    for(int i = 0; i < str.size(); i++){
        // radius[i] = R > i ? min(radius[2*C-i], R-i+1) : 1;
        if(R > i) radius[i] = min(radius[2*C-i], R-i+1);
        else radius[i] = 1;
        while(i + radius[i] < str.size() and i - radius[i] > -1) {
            if(str[i-radius[i]] == str[i+radius[i]] ){
                radius[i] ++;
            }
            else break;
        }
        if(i+radius[i] > R){
            R = i + radius[i]-1;
            C = i;
        }
        Max = max(Max, radius[i]);
    }
    return (Max-1) ;
}
string s;
int main(){
#ifdef DEBUG
    freopen("in", "r", stdin);
#endif  
    int cas = 1;
    while (cin >> s)
    {
        if(s == "END") break;
        printf("Case %d: %d\n", cas++, manacher(s));
    }
    return 0;   
}
```

<h2 id = 12>ac-自动机</h2>


```c++
struct node
{
    int next[26], cnt, fail;
} T[MAX * 10000]; // 要有足够的空间
char s[MAX * 1000];
int tot, rt = 1, cnt;
int Q[MAX*10000], head, tail;
void insert(){
    int tmp = rt;
    for(int i = 0; s[i] != 0; i++){
        if(!T[tmp].next[s[i]-'a']){
            T[tmp].next[s[i]-'a'] = ++tot;
            for(int j = 0; j < 26; j++){
                T[tot].next[j] = 0;
            }
            T[tot].fail = T[tot].cnt = 0;
        }
        tmp = T[tmp].next[s[i]-'a'];
    }
    T[tmp].cnt++;
}
void build(){
    // queue<int> Q;
    // Q.push(rt);
    Q[tail++] = rt;
    while (head < tail)
    {
        int tmp = Q[head++];
        for(int i = 0; i < 26; i++){
            if(T[tmp].next[i]){
                if(tmp == rt){
                    T[T[tmp].next[i]].fail = rt;
                }
                else{
                    int p = T[tmp].fail;
                    while(p){
                        if(T[p].next[i]){
                            T[T[tmp].next[i]].fail = T[p].next[i];
                            break;
                        }
                        p = T[p].fail;
                    }
                    if(!p) T[T[tmp].next[i]].fail = rt;
                }
                // Q.push(T[tmp].next[i]);
                Q[tail++] = T[tmp].next[i];
            }
        }
    }
}
void ac_auto(){
    int p = rt;
    for(int i = 0; s[i] != 0; i++){
        int x = s[i] - 'a';
        while(!T[p].next[x] && p != rt)   
            p = T[p].fail;
        p = T[p].next[x];
        if(!p) p = rt;
        int tmp = p;
        while(tmp != rt){
            if(T[tmp].cnt >= 0){
                cnt += T[tmp].cnt;
                T[tmp].cnt = -1;
            }
            else break;
            tmp = T[tmp].fail;
        }
    }
}
int main()
{
#ifdef DEBUG
    freopen("in.txt", "r", stdin);
#endif
    int cas;
    scanf("%d", &cas);
    while (cas--)
    {
        int n;
        tot = 1;
        head = tail = 0;
        for(int i = 0; i < 26; i++) T[rt].next[i] = 0;
        T[rt].cnt = T[rt].fail = 0;
        cnt = 0;
        scanf("%d", &n);
        while (n--)
        {
            scanf("%s", s);
            insert();
        }
        scanf("%s", s);
        build();
        ac_auto();
        printf("%d\n", cnt);
    }
    return 0;
}
```

Fail 指针和next优化 

> next数组因该改指针数组，提高可读性

```c++
struct node
{
    int next[26], cnt, fail, id;//match;
};

struct AC_auto
{
    int tot, rt;
    node T[MAX];
    void init(){
        tot = 1;
        rt = 1;
        T[rt].cnt = T[rt].fail = 0;
        mset(T[rt].next, 0);
    }
    int insert(char *s, int id){
        int u = rt;
        for(int i = 0; s[i] != 0; i++){
            int c = s[i] - 'a';
            if(!T[u].next[c]){
                T[u].next[c] = ++ tot;
                mset(T[tot].next, 0);
                T[tot].fail = T[tot].cnt = 0;
            } 
            u = T[u].next[c];
        }
        T[u].cnt++;
        T[u].id = id;
        return u;
    }
    void build(){
        queue<int> Q;
        Q.push(rt);
        while(!Q.empty()){
            int u = Q.front();Q.pop();
            for(int i = 0; i < 26; i++){
                if(T[u].next[i]){
                    if(u == rt){
                        T[T[u].next[i]].match = rt;
                        T[T[u].next[i]].fail = rt;
                    }
                    else{
                        int p = T[u].fail;
                        T[T[u].next[i]].fail = T[T[u].fail].next[i];
                        // T[T[u].next[i]].match = T[T[T[u].next[i]].fail].id ? T[T[u].next[i]].fail : T[T[T[u].next[i]].fail].match;
                    }
                    Q.push(T[u].next[i]);
                }
                else{
                    T[u].next[i] = u == rt ? rt : T[T[u].fail].next[i];
                }
            }
        }
    }
    node& operator [] (int idx){
        return T[idx];
    }
}T;

```



<h2 id = 16>回文树</h2>   

```c++
const int MAXN = 10000005;
const int N = 26;
struct Palindromic_Tree
{
	int next[MAXN][N]; //next指针，next指针和字典树类似，指向的串为当前串两端加上同一个字符构成
	int fail[MAXN];	//fail指针，失配后跳转到fail指针指向的节点
	int cnt[MAXN];
	int num[MAXN];
	int len[MAXN]; //len[i]表示节点i表示的回文串的长度
	int S[MAXN];   //存放添加的字符
	int last;	  //指向上一个字符所在的节点，方便下一次add
	int n;		   //字符数组指针
	int p;		   //节点指针

	int newnode(int l)
	{ //新建节点
		for (int i = 0; i < N; ++i)
			next[p][i] = 0;
		cnt[p] = 0;
		num[p] = 0;
		len[p] = l;
		return p++;
	}

	void init()
	{ //初始化
		p = 0;
		newnode(0);
		newnode(-1);
		last = 0;
		n = 0;
		S[n] = -1; //开头放一个字符集中没有的字符，减少特判
		fail[0] = 1;
	}

	int get_fail(int x)
	{ //和KMP一样，失配后找一个尽量最长的
		while (S[n - len[x] - 1] != S[n])
			x = fail[x];
		return x;
	}

	void add(int c)
	{
		c -= 'a';
		S[++n] = c;
		int cur = get_fail(last); //通过上一个回文串找这个回文串的匹配位置
		if (!next[cur][c])
		{											  //如果这个回文串没有出现过，说明出现了一个新的本质不同的回文串
			int now = newnode(len[cur] + 2);		  //新建节点
			fail[now] = next[get_fail(fail[cur])][c]; //和AC自动机一样建立fail指针，以便失配后跳转
			next[cur][c] = now;
			num[now] = num[fail[now]] + 1;
		}
		last = next[cur][c];
		cnt[last]++;
	}

	void count() // 记得跑一遍, 更新cnt
	{
		for (int i = p - 1; i >= 0; --i)
			cnt[fail[i]] += cnt[i];
		//父亲累加儿子的cnt，因为如果fail[v]=u，则u一定是v的子回文串！
	}
};
Palindromic_Tree T;

int main(){
	string s;
	cin >> s;
	T.init();
	for(int i = 0; i < s.size(); i++){
		T.add(s[i]);
	}
	T.count();
	return 0;
}
```
处理从两端插入字符（hdu5421）：维护两个last指针，分为维护前缀和后缀，特别的，如果加入一个字符后整个串是一个回文，那么两个last指针应该只想同一个节点。

<h2 id = 17> 后缀数组 </h2>

```c++
int n, k;
char a[MAX];
int sa[MAX], lcp[MAX], Rank[MAX], tmp[MAX];
bool cmp(int i, int j){
    if(Rank[i] != Rank[j]) return Rank[i] < Rank[j];
    int ri = i + k <= n ? Rank[i+k] : -1;
    int rj = j + k <= n ? Rank[j+k] : -1;
    return ri < rj;
}
int X[MAX], Y[MAX];
void calcu_sa(){
    for(int i = 0; i <= n; i++) {
        sa[i] = i;
        Rank[i] = i < n ? a[i] : -1;
    }
    int m = 200;
    for(int i = 0; i <= m; i++) X[i] = 0;
    for(int i = 0; i <= n; i++) X[Rank[i]+1]++;
    for(int i = 1; i <= m; i++) X[i] += X[i-1];
    for(int i = n; i >= 0; i--) sa[--X[Rank[i]+1]] = i;
    for(k = 1; k <= n; k <<= 1){
        int pos = 0;
        for(int i = n-k+1; i <= n; i++) Y[pos++] = i;
        for(int i = 0; i <= n; i++) if(sa[i] >= k) Y[pos++] = sa[i]-k;
        for(int i = 0; i <= m; i++ ) X[i] = 0;
        for(int i = 0; i <= n; i++) X[Rank[Y[i]]+1]++;
        for(int i = 1; i <= m; i++) X[i] += X[i-1];
        for(int i = n; i >= 0; i--) 
            sa[--X[Rank[Y[i]]+1]] = Y[i];
        tmp[sa[0]] = 0;
        for(int i = 1; i <= n; i++) 
            tmp[sa[i]] = tmp[sa[i-1]] + cmp(sa[i-1], sa[i]);
        for(int i = 0; i <= n; i++) 
            Rank[i] = tmp[i];
        m = Rank[sa[n]]+1;
    }
}
void calcu_lcp(){
    int h = 0;
    lcp[0] = 0;
    for(int i = 0; i < n; i++){
        int j = sa[Rank[i]-1];
        if(h) h--;
        for(; i + h < n and j + h < n; h++)
            if(a[i+h] != a[j+h]) break;
        lcp[Rank[i]-1] = h;
    }
}
```
## SAM

```c++
struct node{
	int len, link, next[26];
};
struct sam{
	node st[MAX];
	int tot, last;
	void init(){
		tot = last = 0;
		st[0].len = 0;
		st[0].link = -1;
		mset(st[0].next, 0);
		tot++;
	}
	void extend(int x){
		int cur = tot++;
		st[cur].len = st[last].len+1;
		int p;
		for(p = last; p!=-1 and !st[p].next[x]; p = st[p].link)
			st[p].next[x] = cur;
		if(p == -1)
			st[cur].link = 0;
		else{
			int q = st[p].next[x];
			if(st[q].len == st[p].len + 1)
				st[cur].link = q;
			else{
				int clone = tot++;
				st[clone].len = st[p].len + 1;
				st[clone].link = st[q].link;
				// memcpy(st[clone].next, st[q].next, sizeof(st[clone].next));
				for(int i = 0; i < 26; i++)
					st[clone].next[i] = st[q].next[i];
				for(; p != -1 and st[p].next[x] == q; p = st[p].link)
					st[p].next[x] = clone;
				st[cur].link = st[q].link = clone;	
			}	
		}
		last = cur;
	}
	int lcs(char *T){
		int u = 0, l = 0, res = 0;
		int n = strlen(T);
		for(int i = 0; i < n; i++){
			int c = T[i]-'a';
			while(u and !st[u].next[c]) 
				u = st[u].link, l = st[u].len;
			if(st[u].next[c])
				u = st[u].next[c], l++;
			res = max(res, l); 
		}
		return res;
	}
}sam;
char s[MAX], t[MAX];
int main(){
#ifdef DEBUG
	freopen64("in", "r", stdin);
#endif
	scanf("%s%s", s, t);
	sam.init();
	for(int i = 0; s[i]; i++)	sam.extend(s[i]-'a');
	printf("%d\n", sam.lcs(t));	
	return 0;
}
```

## GSAM

```c++
struct GSAM
{
    int len[MAX<<1];
    int next[MAX<<1][26];
    int link[MAX<<1];
    int tot;
    void init(){
        tot = 1;
        len[0] = 0;
        link[0] = -1;
    }
    void insert_trie(string &s){
        int rt = 0;
        for(auto c : s) {
            if(next[rt][c-'a']) rt = next[rt][c-'a'];
            else{
                next[rt][c-'a'] = tot++;
                rt = next[rt][c-'a'];
            }
        }
    }
    int insert_sam(int last, int c){
        int cur = next[last][c];
        if(len[cur]) return cur;
        len[cur] = len[last] + 1;
        int p = link[last];
        for(; p != -1 and !next[p][c]; p = link[p]){
            next[p][c] = cur;
        }
        if(p == -1) {
            link[cur] = 0;
            return cur;
        }
        int q = next[p][c];
        if(len[q] == len[p] + 1) {
            link[cur] = q;
            return cur;
        }
        int clone = tot++;
        for(int i = 0; i < 26; i++){
            next[clone][i] = len[next[q][i]] != 0 ? next[q][i] : 0;
        }
        len[clone] = len[p] + 1;
        for(; p != -1 and next[p][c] == q; p = link[p])
            next[p][c] = clone;
        link[clone] = link[q];
        link[cur] = clone;
        link[q] = clone;
        return cur;
    }
    void build(){
        queue<P> Q;
        for(int i = 0; i < 26; i++){
            if(next[0][i]) Q.emplace(i, 0);
        }
        while(!Q.empty()){
            auto cur = Q.front();
            Q.pop();
            auto last = insert_sam(cur.second, cur.first);
            for(int i = 0; i < 26; i++) {
                if(next[last][i]) Q.emplace(i, last);
            }
        }
    }
}sam;

int main(){
#ifdef DEBUG
    freopen64("in", "r", stdin);
#endif
    int T;
    scanf("%d", &T);
    sam.init();
    while (T--)
    {
        string s;
        cin >> s;
        sam.insert_trie(s);
    }
    sam.build();
    LL ans= 0;
    for(int i = 1; i < sam.tot; i++){
        ans += (LL) sam.len[i] - sam.len[sam.link[i]];
    }
    printf("%lld\n", ans);
    return 0;
}
```

---------

<h2 id = 15> FFT </h2>

```c++
int rev[MAX];
void fft(complex<double> *a, int n, int inv){
    int bit = 0; // n==1<<bit
    while((1 << bit) < n) bit++;
    for(int i = 0; i < n; i++){
        rev[i] = (rev[i>>1]>>1) | ((i & 1) << (bit-1));
        if(i < rev[i]) 
            swap(a[i], a[rev[i]]);
    }
    for(int mid = 1; mid < n; mid <<= 1){
        complex<double> wn(cos(PI/mid), inv*sin(PI / mid));
        for(int i = 0; i < n; i += (mid<<1)){
            complex<double> w(1, 0);
            for(int j = 0; j < mid; j++, w *= wn){
                complex<double> u = a[i+j], t = w * a[i+j+mid];
                a[i+j] = u + t, a[i+j+mid] = u-t;
            }
        }
    }
}
complex<double> a[MAX], b[MAX];
int main(){
#ifdef DEBUG
    freopen("in.txt", "r", stdin);
#endif
    int n, m;
    cin >> n >> m;
    for(int i = 0; i <= n; i++){
        double x;
        scanf("%lf", &x);
        a[i] = complex<double>(x, 0);
    }
    for(int i = 0; i <= m; i++) {
        double x;
        scanf("%lf", &x);
        b[i] = complex<double>(x, 0);
    }
    int N = 1;
    while(N <= n + m) N <<= 1;
    fft(a, N, 1); fft(b, N, 1);
    for(int i = 0; i <= N; i++){
        a[i] *= b[i];
    }
    fft(a, N, -1);
    for(int i = 0; i <= n+m; i++){
        printf("%d ", (int)(a[i].real() / N + 0.5));
    }
    return 0;
}
```

----------
<h2 id = 3> 快读，快写 </h2>  

```c++
inline ll read()
{
    ll ret = 0, flag = 1;
    char c = getchar();
    while (c < '0' || c > '9')
    {
        if (c == '-')
            flag = -1;
        c = getchar();
    }
    while (c >= '0' && c <= '9')
        ret = 10 * ret + (c ^ '0'), c = getchar();
    return flag * ret;
}
inline void write(ll x)
{
    if (x < 0)
        putchar('-'), x = -x;
    if (x > 9)
        write(x / 10);
    putchar(x % 10 + '0');
}
```


// Rigorous evaluation of the power sums in the tile-moment condition (C), base 3/2.
// Integer interval arithmetic at scale 2^42 with outward rounding (primitives as in ../exact/check.cpp).
//
// For each mesh interval J and each upper-path prefix p, over all 3^T residue tiles t and all
// phase pieces of J (maximal subintervals on which the T phase maps are affine), computes an
// upper bound A_{t,piece} for  sum_{a<T} wgt[a] * sup_{x in piece} |G_a(x)| , where G_a is the
// complex sum of the paths from the top of p that end at age a, and wgt[a] >= q^T y_a c(y_a).
// Remove the common critical suffix of each age group before evaluating its modulus.
// Bound the resulting return sum by its endpoint moduli and the chord error
// h^2 sup|H''|/8, where prep_tile.py provides rational curvature bounds.
// Prints upper bounds for P_j = sum A^j, j=0..M,
// and E(J) = sum over tiles and pieces of the multiplicity mu = number of cells met by the
// image of the piece (the induction needs the sup over that union).  P_j is weighted by mu.
//
// Input (stdin, produced by prep_tile.py):
//   T M KX nthreads
//   twopi_lo twopi_hi
//   wgt_0 ... wgt_{T-1}
//   critical digits d_1 ... d_T
//   NB  then NB cell boundaries X (cells [X_k,X_{k+1})/2^KX partition [-1/2,1/2))
//   NP  then NP lines  prefix-digit-string (length T), then NP lines of T curvature bounds
// Output: lines "E j E(J)" and "P j p P_0 ... P_M" (P_j scaled by 2^42).
#include <algorithm>
#include <atomic>
#include <bit>
#include <cassert>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <limits>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#ifdef CERT_INDEPENDENT_ARITHMETIC
#include "../independent_interval.hpp"
#endif

using Int = std::int64_t;
using Wide = __int128_t;
constexpr Int scale = Int{1} << 42;

void require(bool condition, const char* message) {
  if (!condition) {
    std::fprintf(stderr, "tilecheck: %s\n", message);
    std::abort();
  }
}
Int narrow(Wide x){
  require(x >= std::numeric_limits<Int>::min() && x <= std::numeric_limits<Int>::max(),
          "signed 64-bit range exceeded");
  return (Int)x;
}
Int floor_div(Wide x, Int dv){ assert(dv>0); return narrow(x/dv - (x%dv<0)); }
Int ceil_div(Wide x, Int dv){ assert(dv>0); return narrow(x/dv + (x%dv>0)); }
Wide floor_divw(Wide x, Wide dv){ return x/dv - (x%dv<0); }
struct Interval{ Int lo=0, hi=0; };
Interval operator+(Interval a, Interval b){ return {narrow(Wide(a.lo)+b.lo), narrow(Wide(a.hi)+b.hi)}; }
Interval operator-(Interval a){ return {narrow(-Wide(a.hi)), narrow(-Wide(a.lo))}; }
Interval operator*(Interval a, Interval b){
  Wide p[4]={Wide(a.lo)*b.lo, Wide(a.lo)*b.hi, Wide(a.hi)*b.lo, Wide(a.hi)*b.hi};
  return {floor_div(*std::min_element(p,p+4),scale), ceil_div(*std::max_element(p,p+4),scale)};
}
Int mul_up(Int a, Int b){ assert(a>=0&&b>=0); return ceil_div(Wide(a)*b,scale); }
struct Complex{ Interval re, im; };
Complex operator*(Complex a, Complex b){ return {a.re*b.re + -(a.im*b.im), a.re*b.im + a.im*b.re}; }
Complex operator+(Complex a, Complex b){ return {a.re+b.re, a.im+b.im}; }
bool isZero(const Complex& z){ return z.re.lo==0&&z.re.hi==0&&z.im.lo==0&&z.im.hi==0; }
Complex exp_unit(Interval angle){   // |angle| <= pi < 22/7; prep checks (22/7)^32/32! < 2^-42
#ifdef CERT_INDEPENDENT_ARITHMETIC
  const auto value = crosscheck::horner_exponential(crosscheck::widen_input(angle.lo,angle.hi));
  const auto re = crosscheck::shrink_output(value.re);
  const auto im = crosscheck::shrink_output(value.im);
  return {{re.lo,re.hi},{im.lo,im.hi}};
#else
  Complex out{{scale,scale},{0,0}}; Interval term{scale,scale};
  for(Int n=1;n<32;++n){ term=term*angle; term={floor_div(term.lo,n),ceil_div(term.hi,n)};
    switch(n%4){ case 0: out.re=out.re+term; break; case 1: out.im=out.im+term; break;
                 case 2: out.re=out.re+-term; break; case 3: out.im=out.im+-term; break; } }
  return {{narrow(Wide(out.re.lo)-1),narrow(Wide(out.re.hi)+1)},
          {narrow(Wide(out.im.lo)-1),narrow(Wide(out.im.hi)+1)}};
#endif
}
Int magnitude_up(const Complex& z){
#ifdef CERT_INDEPENDENT_ARITHMETIC
  return crosscheck::modulus_up({{z.re.lo,z.re.hi},{z.im.lo,z.im.hi}});
#else
  const auto abs_wide = [](Int v) -> Wide { return v < 0 ? -Wide(v) : Wide(v); };
  const Wide re=std::max(abs_wide(z.re.lo),abs_wide(z.re.hi));
  const Wide im=std::max(abs_wide(z.im.lo),abs_wide(z.im.hi));
  require(re <= (Wide(1)<<60) && im <= (Wide(1)<<60), "modulus square range exceeded");
  const Wide sq=Wide(re)*re+Wide(im)*im; if(sq==0) return 0;
  const auto high=(std::uint64_t)(sq>>64);
  const int bits= high? 64+std::bit_width(high) : std::bit_width((std::uint64_t)sq);
  assert(bits<=122); Int root=Int{1}<<((bits+1)/2);
  while(true){ const Int nx=(Int)((root+sq/root)/2); if(nx>=root) return root+(Wide(root)*root!=sq); root=nx; }
#endif
}

int T,M,KX,NT; Int twopi_lo,twopi_hi; std::vector<Int> wgt; std::vector<int> dig;
std::vector<std::string> pre; int NP; Wide DEN; Wide TH;   // DEN=2^KX 3^T, TH=DEN/4
std::vector<Wide> slope;                                 // slope of Phi_i in X: 2^i 3^{T-i}
std::vector<std::vector<Int>> curv;                      // per prefix, per age: scaled curvature bound
const Complex ONE{{scale,scale},{0,0}}, ZERO{{0,0},{0,0}};

Complex expPhase(Wide Phi){       // e(Phi/DEN)
  require(-DEN/2 <= Phi && Phi <= DEN/2, "phase outside centered interval");
  Wide a = (Phi>=0)? Wide(twopi_lo)*Phi : Wide(twopi_hi)*Phi;
  Wide b = (Phi>=0)? Wide(twopi_hi)*Phi : Wide(twopi_lo)*Phi;
  Interval ang{narrow(floor_divw(a,DEN)), narrow(-floor_divw(-b,DEN))};
  return exp_unit(ang);
}
struct End {
  std::vector<Wide> Phi;
  std::vector<Complex> twists;
  std::vector<Complex> returns, top;
  std::vector<Int> return_norms;

  End() : returns(NP*T, ZERO), top(NP, ONE), return_norms(NP*T, 0) {}
};
struct Piece{ Wide Xa,Xb; End A,B; };

// Compute the sum of paths returning at this step. If the preceding
// return was r steps earlier, the intervening critical prefix is fixed.
void stepAll(End& endpoint, int step, const Complex& twist) {
  std::vector<std::pair<int, Complex>> excursions;
  for (int r=2; r<step; ++r) {
    if (dig[r-1] == 0) continue;
    Complex factor = ONE;
    for (int j=0; j<r-1; ++j) {
      if (dig[j] == 1) factor = factor * endpoint.twists[step-r+j];
    }
    excursions.emplace_back(step-r-1, factor);
  }
  for (int p=0; p<NP; ++p) {
    Complex* history = &endpoint.returns[p*T];
    Complex total = step == 1 ? ZERO : history[step-2]; // return immediately by digit 0
    for (const auto& [previous, factor] : excursions) {
      if (!isZero(history[previous])) total = total + history[previous] * factor;
    }
    if (pre[p][step-1] == '1') {
      total = total + endpoint.top[p];
      endpoint.top[p] = endpoint.top[p] * twist;
    }
    history[step-1] = total;
    endpoint.return_norms[p*T+step-1] = magnitude_up(total);
  }
  endpoint.twists.push_back(twist);
}
// rebuild an endpoint from scratch given its phase history Phi[0..i]
End rebuild(const std::vector<Wide>& Phi){
  End E; E.Phi=Phi;
  for(size_t s=1;s<Phi.size();s++) stepAll(E,(int)s,expPhase(Phi[s]));
  return E;
}
struct Acc{ std::vector<std::vector<Wide>> P; long long E=0; };

std::vector<Wide> cellB;   // cell boundaries times 3^T (comparable with Phi = phi*DEN)
int mult(Wide u, Wide v){  // #cells [c_k,c_{k+1}) meeting [u,v)
  int k0=int(std::upper_bound(cellB.begin(),cellB.end(),u)-cellB.begin())-1;   // cell containing u
  int k1=int(std::lower_bound(cellB.begin(),cellB.end(),v)-cellB.begin())-1;   // last cell with c_k < v
  require(u<v && k0>=0 && k1>=k0 && k1<(int)cellB.size()-1,
          "invalid image interval");
  return k1-k0+1;
}
void leaf(const std::vector<Piece>& pcs, Acc& acc){
  for(const Piece& pc:pcs){
    const int mu=mult(pc.A.Phi[T],pc.B.Phi[T]);
    require(!__builtin_add_overflow(acc.E, mu, &acc.E), "exit-cell count overflow");
    Wide h=pc.Xb-pc.Xa;  // chord factor h^2/8 in units of 2^{-2KX}
    for(int p=0;p<NP;p++){
      Int A=0;
      for(int a=0;a<T;a++){
        // A path ending at age a is a return at T-a followed by a
        // common critical suffix of unit modulus.
        const int index = p*T + T-a-1;
        Int m=std::max(pc.A.return_norms[index],pc.B.return_norms[index]);
        // err = curv * h^2 / (8 * 2^{2KX})
        Int err=narrow(-floor_divw(-(Wide(curv[p][a])*h*h), Wide(8)<<(2*KX)));
        A=narrow(Wide(A)+mul_up(wgt[a],narrow(Wide(m)+err)));
      }
      Int pw=scale;
      for(int j=0;j<=M;j++){
        auto& sum = acc.P[p][j];
        require(!__builtin_add_overflow(sum, Wide(mu)*pw, &sum), "power sum overflow");
        if(j<M) pw=mul_up(pw,A);
      }
    }
  }
}
void dfs(int i, const std::vector<Piece>& pcs, Acc& acc){
  if(i==T){ leaf(pcs,acc); return; }
  for(int r=0;r<3;r++){
    std::vector<Piece> nx;
    for(const Piece& pc:pcs){
      Wide a=pc.A.Phi[i], b=pc.B.Phi[i];    // phi_i at the ends (increasing in X)
      auto emit=[&](Wide Xa,Wide Xb,const End& EA,const End& EB,int c){
        Piece q{Xa, Xb, EA, EB};
        for(End* E:{&q.A,&q.B}){
          Wide ph=E->Phi[i]; Wide num=2*ph+Wide(c)*DEN;
          require(num%3==0, "phase endpoint not on exact denominator");
          E->Phi.push_back(num/3);
        }
        Complex eA=expPhase(q.A.Phi[i+1]), eB=expPhase(q.B.Phi[i+1]);
        stepAll(q.A,i+1,eA); stepAll(q.B,i+1,eB); nx.push_back(std::move(q)); };
      auto cutAt=[&](Wide target)->Wide{
        Wide num=target-a;
        require(num%slope[i]==0, "cut not on exact dyadic grid");
        const Wide cut=pc.Xa+num/slope[i];
        require(pc.Xa<cut && cut<pc.Xb, "cut outside the interior of its piece");
        return cut;
      };
      auto midEnd=[&](Wide Xc){ std::vector<Wide> Ph(i+1); for(int s=0;s<=i;s++) Ph[s]=pc.A.Phi[s]+slope[s]*(Xc-pc.Xa); return rebuild(Ph); };
      if(r==0){ emit(pc.Xa,pc.Xb,pc.A,pc.B,0); continue; }
      if(r==1){ // c=1 if phi<1/4 else -2
        if(b<=TH) emit(pc.Xa,pc.Xb,pc.A,pc.B,1);
        else if(a>=TH) emit(pc.Xa,pc.Xb,pc.A,pc.B,-2);
        else { Wide Xc=cutAt(TH); End Ec=midEnd(Xc); emit(pc.Xa,Xc,pc.A,Ec,1); emit(Xc,pc.Xb,Ec,pc.B,-2); }
      } else { // r==2: c=-1 if phi>=-1/4 else 2
        if(a>=-TH) emit(pc.Xa,pc.Xb,pc.A,pc.B,-1);
        else if(b<=-TH) emit(pc.Xa,pc.Xb,pc.A,pc.B,2);
        else { Wide Xc=cutAt(-TH); End Ec=midEnd(Xc); emit(pc.Xa,Xc,pc.A,Ec,2); emit(Xc,pc.Xb,Ec,pc.B,-1); }
      }
    }
    dfs(i+1,nx,acc);
  }
}
std::string w128(Wide x){ require(x>=0,"negative power sum"); if(x==0) return "0"; std::string s; while(x>0){ s+=char('0'+int(x%10)); x/=10; } return std::string(s.rbegin(),s.rend()); }

int main(int argc, char** argv){
  const bool no_symmetry = argc == 2 && std::string(argv[1]) == "--no-symmetry";
  if (argc != 1 && !no_symmetry) {
    std::fprintf(stderr, "usage: tilecheck [--no-symmetry]\n");
    return 2;
  }
  std::cin>>T>>M>>KX>>NT>>twopi_lo>>twopi_hi;
  require(std::cin && 1<=T && T<=20 && 1<=M && M<=64 && T+2<=KX && KX<=24
          && 1<=NT && NT<=256 && 0<twopi_lo && twopi_lo<=twopi_hi && twopi_hi<8*scale,
          "unsupported checker dimensions or pi bounds");
  wgt.resize(T); for(auto&w:wgt) { std::cin>>w; require(std::cin && w>=0,"invalid weight"); }
  dig.resize(T); for(auto&x:dig) { std::cin>>x; require(std::cin && (x==0 || x==1),"invalid digit"); }
  require(dig[0]==1,"invalid first critical digit");
  int NB=0; std::cin>>NB;
  require(std::cin && 2<=NB && NB<=65536,"invalid boundary count");
  std::vector<Int> bnd(NB); for(auto&b:bnd) std::cin>>b;
  const Int half_grid = Int(1)<<(KX-1);
  require(std::cin && bnd.front()==-half_grid && bnd.back()==half_grid
          && std::adjacent_find(bnd.begin(),bnd.end(),std::greater_equal<Int>())==bnd.end(),
          "phase cells must partition the centered interval");
  int NJ=NB-1; std::vector<std::pair<Int,Int>> Js(NJ); for(int j=0;j<NJ;j++) Js[j]={bnd[j],bnd[j+1]};
  bool symmetric = !no_symmetry;
  for (int i=0; i<NB; ++i) symmetric &= Wide(bnd[i]) == -Wide(bnd[NB-1-i]);
  // Reflection negates the residue word and conjugates the return sums
  // on open pieces. Their suprema and the number of image cells agree.
  const int jobs = symmetric ? (NJ+1)/2 : NJ;
  std::cin>>NP;
  require(std::cin && 1<=NP && NP<=(1<<20),"invalid prefix count");
  pre.resize(NP);
  for(auto&s:pre){
    std::cin>>s;
    require(std::cin && s.size()==static_cast<size_t>(T) && s.find_first_not_of("01")==std::string::npos,
            "invalid upper prefix");
  }
  Wide p3=1; for(int i=0;i<T;i++) p3*=3; DEN=(Wide(1)<<KX)*p3; TH=DEN/4; assert(TH*4==DEN);
  for(Int b:bnd) cellB.push_back(Wide(b)*p3);
  slope.resize(T+1); for(int i=0;i<=T;i++){ Wide s=Wide(1)<<i; for(int k=0;k<T-i;k++) s*=3; slope[i]=s; }
  // curvature bounds per prefix and age: ceil(2^42 (2pi)^2 sum_paths theta'^2), from prep_tile.py
  curv.assign(NP,std::vector<Int>(T,0));
  for(int p=0;p<NP;p++) for(int a=0;a<T;a++){
    std::cin>>curv[p][a];
    require(std::cin && curv[p][a]>=0,"invalid curvature bound");
  }
  std::vector<Acc> res(jobs); std::atomic<int> next{0}; std::mutex mu;
  auto worker=[&](){ while(true){ int j=next++; if(j>=jobs) break;
      Acc acc; acc.P.assign(NP,std::vector<Wide>(M+1,0));
      Piece pc; pc.Xa=Js[j].first; pc.Xb=Js[j].second;
      pc.A.Phi={Wide(pc.Xa)*p3}; pc.B.Phi={Wide(pc.Xb)*p3};
      dfs(0,{pc},acc); res[j]=std::move(acc);
      std::lock_guard<std::mutex> g(mu); std::fprintf(stderr,"J %d/%d done\n",j+1,jobs); } };
  std::vector<std::thread> th; for(int k=0;k<NT;k++) th.emplace_back(worker); for(auto&t:th) t.join();
  for(int j=0;j<NJ;j++){
    const auto& acc = res[symmetric ? std::min(j, NJ-1-j) : j];
    std::printf("E %d %lld\n",j,acc.E);
    for(int p=0;p<NP;p++){ std::printf("P %d %d",j,p); for(int k=0;k<=M;k++) std::printf(" %s",w128(acc.P[p][k]).c_str()); std::printf("\n"); } }
}

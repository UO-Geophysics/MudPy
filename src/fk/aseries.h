c global parameters shared by routines
      INTEGER max_ray_seg,max_num_lay,max_num_pts,max_terms,mfft
      REAL    pi
      PARAMETER (max_ray_seg=5000, max_num_lay=2000, mfft=2048,
     &           max_num_pts=5000, max_terms=10, pi=3.141592653)
c	max_ray_seg---- maximum number of ray segments
c	max_num_lay---- maximum number of layers in the model
c	max_num_pts---- maximum number of points in time series
c	max_terms------ maximum number of terms for Kn expansion
c	mfft----------- maximum number of fft when qray


c Block for expansion of Bessel function
      REAL	kfn(5, max_terms)	! Kn(x), n=0,1,2 and nKn(x)/x
      INTEGER	order			! # of expansion terms
      COMMON	/kfnct/order, kfn


c Block for velocity model
      LOGICAL	debug, tele
      INTEGER num_lay, src_lay, src, nCom
      REAL vps(2,max_num_lay),den(max_num_lay),
     +thk(max_num_lay),q(2,max_num_lay)
      REAL Vmin, Hmin
      PARAMETER (Vmin=0.001, Hmin=0.01)
      COMMON/model/num_lay,vps,den,thk,q,src_lay,src,nCom,
     +		   debug,tele
c 	num_lay		number of layers in the model
c	vps		p and s velocities
c	den		density
c	q		Qp and Qs
c	thk		thickness of the layer
c	src_lay		the layer where source is located at the top.
c	src		0=> explosion; 1=> single force; 2=>double-couple
c	nCom		number of components of Green's function (<=9)


c Block for ray description
      CHARACTER*1 src_type, seg_ps(max_ray_seg)
      INTEGER num_ray_seg, bttm, topp, ray_seg(max_ray_seg)
      REAL ray_len(2,max_num_lay), dt20, a0, a1, t0, tc, ts
      COMPLEX p0, pc, coef0(9,max_terms)
      COMMON/ray/num_ray_seg,ray_seg,ray_len,bttm,topp,dt20,a0,a1,
     +		t0,tc,ts,p0,pc,coef0,src_type,seg_ps
c	num_ray_seg	number of ray segments for this ray
c	ray_seg		layer the segment is in
c	seg_ps		p or s wave for this segment
c	topp,bttm	the shallowest/deepest layer the ray has penetrated
c	ray_len		the total vertical length of ray in each layer
c	src_type	source type (p, s, t)
c	t0, p0		at first-motion arrival time
c	tc, pc		at head-wave arrival time
c	ts		t* = T/Q
c	dt20		abs of 2nd and 3rd direvatives of t(p) at p0
c	a0, a1		at t=t0, dp/dt = a0+i*a1/sqrt(t-t0)
c	coef0		refl/trans. coef. (including src and recv) at t0


c Block for response traces
      INTEGER num_pts
      REAL dt,tstart,tend
      REAL ph0(max_num_pts,9,max_terms),ph1(max_num_pts,9,max_terms)
      COMMON/rsl/num_pts,dt,tstart,tend,ph0,ph1
c	num_pts		number of points of time series
c	dt		sampling rate
c	tstart		starting time
c	tend		ending time
c	ph0		components of green's function (first motion appr.)
c	ph1		components of green's function (first+remaining term)

c***************************************************************
c  F-K: @(#) model.h			1.0 4/29/2000
c
c  Copyright (c) 2000 by L. Zhu
c  See README file for copying and redistribution conditions.
c
c Inlucde header file pertaining the velocity model info
c       velocity model in common block /model/:
c               mb---  number of layers.
c               src--- source layer index (on top of it).
c		rcv--- receiver layer index (on top of it).
c		stype- source type (0=ex;1=sf;2=dc).
c               ka--- (w/vp)^2.
c               kb--- (w/vs)^2.
c               d---- thickness of the layer.
c               rho-- density of the layer.
c		mu--- shear modulus.
c		xi--- mu/bulk_modulus.
c		si(3,6) ---- source coefs. of n=0,1,2.
      integer mb,stype,src,rcv,nlay,ndis,nt,updn
c max. # of layers and receivers and time length
      parameter(nlay=20, ndis=5000, nt=4096)
      real d(nlay),rho(nlay),mu(nlay),xi(nlay),si(3,6),epsilon
      complex ka(nlay),kb(nlay)
      PARAMETER (epsilon=0.0001)
      common /model/mb,stype,src,rcv,updn,ka,kb,d,rho,mu,xi,si

*
*  -------------------------------------------------------
*
*       qgraf-3.3
*       a module for Feynman diagram generation
*
*       copyright 1990-2018 by P. Nogueira
*
*       reference:
*        [1] Automatic Feynman graph generation
*            P. Nogueira
*            J. Comput. Phys. 105 (1993) 279-289
*            https://doi.org/10.1006/jcph.1993.1074
*
*       documentation:
*         [2] files 'qgraf-3.0.pdf' and 'qgraf-3.3.pdf'
*
*  -------------------------------------------------------
*
      program qgraf
      implicit integer(a-z)
      save
      parameter ( maxleg=11, maxrho=7, maxi=11, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( sxbuff=2040 )
      common/z1in/leg(1:maxleg),nleg,nloop,incom
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1
      common/z2g/dis,dsym
      common/z3g/nv(0:0),dpntro(0:0),vparto(0:0),vdeg(0:0),nvert
      common/z4g/n,nli
      common/z5g/psym(0:0),psyms,nsym
      common/z6g/p1(1:maxleg),invp1(1:maxleg)
      common/z7g/lmap(1:maxn,1:maxdeg),vmap(1:maxn,1:maxdeg),
     :pmap(1:maxn,1:maxdeg),vlis(1:maxn),invlis(1:maxn)
      common/z8g/degree(1:maxn),xn(1:maxn)
      common/z9g/tpc(0:0)
      common/z11g/nphi,nblok,nprop,npprop
      common/z12g/jflag(1:11),mflag(1:24)
      common/z14g/zcho(0:maxli),zbri(0:maxli),zpro(0:maxli),
     :rbri(0:maxli),sbri(0:maxli)
      common/z16g/rdeg(1:maxn),amap(1:maxn,1:maxdeg)
      common/z17g/xtail(1:maxn),xhead(1:maxn),ntadp
      common/z18g/eg(1:maxn,1:maxn),flow(1:maxli,0:maxleg+maxrho)
      common/z19g/vfo(1:maxn)
      common/z20g/tftyp(0:0),tfnarg(0:0),tfa(0:0),tfb(0:0),tfc(0:0),
     :tfo(0:0),tf2(0:0),ntf
      character*(srec) auxlin
      common/z22g/auxlin
      common/z25g/ex(1:maxli),ey(1:maxli),ovm(1:maxn,1:maxdeg)
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z34g/link(0:0),antiq(0:0),blok(0:0)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z43g/vmkr(0:0),vmkmao(0:0),vmkmio(0:0)
      integer xli(2*maxn+2*maxrho+maxdeg)
      call init0
      jflag(6)=1
      i=maxdeg
      if(i.lt.3)then
      auxlin(1:srec)='parameter "maxdeg" is not properly set'
      call messag(1,0,0,0)
      endif
      i=maxleg
      if(i.lt.3)then
      auxlin(1:srec)='parameter "maxleg" is not properly set'
      call messag(1,0,0,0)
      endif
      i=maxrho
      if(i.lt.3)then
      auxlin(1:srec)='parameter "maxrho" is not properly set'
      call messag(1,0,0,0)
      endif
      i=maxi
      j=max(maxleg,maxrho)
      if(i.lt.j)then
      auxlin(1:srec)='parameter "maxi" is not properly set'
      call messag(1,0,0,0)
      endif
      i=maxn
      j=2*max(maxleg,maxrho)-2
      if(i.lt.j)then
      auxlin(1:srec)='parameter "maxn" is not properly set'
      call messag(1,0,0,0)
      endif
      i=maxli
      j=2*max(maxleg,maxrho)+maxrho-3
      if(i.lt.j)then
      auxlin(1:srec)='parameter "maxli" is not properly set'
      call messag(1,0,0,0)
      endif
      call spp(0)
      call inputs
      if(jflag(1).ne.0)then
      call iki
      call prepos(1)
      endif
      call rflag
      if(jflag(2).ne.0)then
      goto 94
      endif
      rho(1)=nleg
      rho(2)=0
      rhop1=rho(1)+1
      jj=rhop1-3+2*nloop
      do i1=nrho,3,-1
      rho(i1)=jj/(i1-2)
      jj=jj-rho(i1)*(i1-2)
      enddo
      goto 50
   18 continue
      do i1=4,nrho
      if(rho(i1).gt.0)then
      goto 30
      endif
      enddo
      goto 94
   30 continue
      rho(i1)=rho(i1)-1
      jj=i1-2+rho(3)
      do i2=i1-1,3,-1
      rho(i2)=jj/(i2-2)
      jj=jj-rho(i2)*(i2-2)
      enddo
   50 continue
      do i1=3,nrho
      if(stib(nv(0)+i1).eq.0)then
      if(rho(i1).gt.0)then
      goto 18
      endif
      endif
      enddo
      jflag(3)=0
      nili=-rho(1)
      do i1=3,nrho
      nili=nili+i1*rho(i1)
      enddo
      nili=nili/2
      call sdiag(1,-1)
      cntr10=0
      if(npprop.eq.0)then
      if(nili.gt.0)then
      goto 41
      endif
      endif
      if(zpro(nili).eq.0)then
      goto 41
      endif
      do i1=nloop,nili
      if(zcho(i1).ne.0)then
      if(zbri(nili-i1).ne.0)then
      do i2=0,nili-i1
      if(rbri(i2).ne.0)then
      if(sbri(nili-i1-i2).ne.0)then
      goto 89
      endif
      endif
      enddo
      endif
      endif
      enddo
      goto 41
   89 continue
      if(mflag(20).ne.0)then
      do i1=1,ntf
      if(abs(stib(tftyp(0)+i1)).eq.6)then
      j1=0
      j2=0
      jj=stib(stib(tfo(0)+i1)+1)
      if(stib(vmks(0)+jj).ne.1)then
      auxlin(1:srec)='main_1'
      call messag(1,0,0,0)
      endif
      do i2=nrho,3,-1
      if(rho(i2).gt.0)then
      j1=j1+rho(i2)*stib(stib(vmkmio(0)+jj)+i2)
      j2=j2+rho(i2)*stib(stib(vmkmao(0)+jj)+i2)
      endif
      enddo
      if(stib(tftyp(0)+i1).gt.0)then
      if(j1.gt.stib(tfb(0)+i1))then
      goto 41
      elseif(j2.lt.stib(tfa(0)+i1))then
      goto 41
      endif
      else
      if(j1.ge.stib(tfa(0)+i1))then
      if(j2.le.stib(tfb(0)+i1))then
      goto 41
      endif
      endif
      endif
      endif
      enddo
      endif
      cntr10=1
   40 continue
      call gen10(cntr10)
      goto 42
   41 continue
      jflag(3)=1
   42 continue
      if(cntr10.eq.0)then
      call spp(2)
      goto 18
      endif
      cntr10=0
      do i1=1,rho(1)
      pmap(i1,1)=leg(p1(i1))
      pmap(vmap(i1,1),lmap(i1,1))=stib(link(0)+leg(p1(i1)))
      enddo
      vind=rho(1)
   57 continue
      vind=vind+1
      vv=vlis(vind)
      vfo(vv)=stib(stib(dpntro(0)+degree(vv))+pmap(vv,1))
  100 continue
      do i1=1,rdeg(vv)
      ii=stib(vfo(vv)+i1)-pmap(vv,i1)
      if(ii.gt.0)then
      goto 104
      elseif(ii.lt.0)then
      vfo(vv)=vfo(vv)+degree(vv)+1
      goto 100
      endif
      enddo
      goto 163
   58 continue
      vfo(vv)=vfo(vv)+degree(vv)+1
      do i1=1,rdeg(vv)
      if(stib(vfo(vv)+i1).ne.pmap(vv,i1))then
      goto 104
      endif
      enddo
      goto 163
  104 continue
      if(vind.ne.rhop1)then
      vind=vind-1
      vv=vlis(vind)
      goto 58
      endif
      goto 40
  163 continue
      sdeg=rdeg(vv)+g(vv,vv)
      do i1=rdeg(vv)+1,sdeg,2
      j1=vfo(vv)+i1
      j2=stib(j1)
      j3=stib(j1+1)
      if(j2.gt.j3)then
      goto 58
      elseif(j2.ne.stib(link(0)+j3))then
      goto 58
      elseif(i1+1.ne.sdeg)then
      if(j2.gt.stib(j1+2))then
      goto 58
      endif
      endif
      enddo
      do i1=rdeg(vv)+1,degree(vv)
      pmap(vv,i1)=stib(vfo(vv)+i1)
      enddo
      do i1=sdeg+1,degree(vv)-1
      if(vmap(vv,i1).eq.vmap(vv,i1+1))then
      if(pmap(vv,i1).gt.pmap(vv,i1+1))then
      goto 58
      endif
      endif
      enddo
      if(mflag(18).ne.0)then
      do i1=rdeg(vv)+1,degree(vv)
      if(stib(tpc(0)+pmap(vv,i1)).eq.5)then
      goto 58
      endif
      enddo
      endif
      if(mflag(19).ne.0)then
      do i1=1,ntadp
      if(xtail(i1).eq.vv)then
      jj=xhead(i1)
      elseif(xhead(i1).eq.vv)then
      jj=xtail(i1)
      else
      jj=0
      endif
      if(jj.ne.0)then
      do i2=1,degree(vv)
      if(jj.eq.vmap(vv,i2))then
      if(stib(tpc(0)+pmap(vv,i2)).eq.1)then
      goto 58
      endif
      endif
      enddo
      endif
      enddo
      endif
      do i1=sdeg+1,degree(vv)
      pmap(vmap(vv,i1),lmap(vv,i1))=stib(link(0)+pmap(vv,i1))
      enddo
      if(vind.lt.n)then
      goto 57
      endif
      if(mflag(12).eq.0)then
      goto 333
      endif
      if(nloop.eq.0)then
      if(mflag(12).gt.0)then
      goto 333
      elseif(mflag(12).lt.0)then
      goto 58
      endif
      endif
      do i1=1,n
      xli(i1)=0
      enddo
      do 74 i1=1,n
      if(xli(i1).eq.0)then
      ii=i1
      kk=1
  820 continue
      do i2=1,degree(ii)
      if(i2.ne.xli(ii))then
      jj=pmap(ii,i2)
      if(stib(antiq(0)+jj).ne.0)then
      k=ii
      ii=vmap(k,i2)
      xli(ii)=lmap(k,i2)
      goto 830
      endif
      endif
      enddo
      goto 74
  830 continue
      if(ii.gt.rho(1))then
      if(ii.ne.i1)then
      kk=-kk
      goto 820
      endif
      if(kk.gt.0)then
      if(mflag(12).gt.0)then
      goto 58
      else
      goto 333
      endif
      endif
      endif
      endif
   74 continue
      if(mflag(12).lt.0)then
      goto 58
      endif
  333 continue
      if(jflag(11).ne.0)then
      dsym=nsym
      else
      dsym=1
      jk=psym(0)-rho(1)
      do 206 i1=2,nsym
      if(mflag(14).eq.0)then
      do i2=rhop1,n
      j1=stib(jk+i2)
      j2=rdeg(i2)+1
  400 continue
      if(j2.le.degree(i2))then
      j3=vmap(i2,j2)
      j4=1
  410 continue
      if(j4.le.degree(i2))then
      if(vmap(j1,j4).ne.stib(jk+j3))then
      j4=j4+1
      goto 410
      endif
      endif
      do i3=1,g(i2,j3)
      xli(i3)=pmap(j1,j4+i3-1)
      enddo
      do i3=1,g(i2,j3)-1
      do i4=i3+1,g(i2,j3)
      if(xli(i3).gt.xli(i4))then
      ii=xli(i3)
      xli(i3)=xli(i4)
      xli(i4)=ii
      endif
      enddo
      enddo
      do i3=1,g(i2,j3)
      ii=xli(i3)-pmap(i2,j2)
      if(ii.lt.0)then
      goto 58
      elseif(ii.gt.0)then
      goto 339
      endif
      j2=j2+1
      enddo
      goto 400
      endif
      enddo
      endif
      if(mflag(17).ne.0)then
      do i2=rhop1,n
      j1=stib(jk+i2)
      if(i2.ne.j1)then
      ii=stib(vfo(j1))-stib(vfo(i2))
      if(ii.lt.0)then
      goto 58
      elseif(ii.gt.0)then
      goto 339
      endif
      endif
      enddo
      endif
      dsym=dsym+1
  339 continue
      jk=jk+(n-rho(1))
  206 continue
      endif
      do 314 i=1,ntf
      if(stib(tfnarg(0)+i).eq.0)then
      auxlin(1:srec)='main_2'
      call messag(1,0,0,0)
      endif
      if(abs(stib(tftyp(0)+i)).gt.10)then
      goto 314
      endif
      if(abs(stib(tftyp(0)+i)).eq.1)then
      ii=0
      do 319 j=rhop1,n
      do jj=rdeg(j)+1,rdeg(j)+g(j,j),2
      do ij=1,stib(tfnarg(0)+i)
      if(stib(stib(tfo(0)+i)+ij).eq.pmap(j,jj))then
      ii=ii+1
      elseif(stib(link(0)+stib(stib(tfo(0)+i)+ij)).eq.
     :pmap(j,jj))then
      ii=ii+1
      endif
      enddo
      enddo
      do jj=rdeg(j)+g(j,j)+1,degree(j)
      do ij=1,stib(tfnarg(0)+i)
      if(stib(stib(tfo(0)+i)+ij).eq.pmap(j,jj))then
      ii=ii+1
      elseif(stib(link(0)+stib(stib(tfo(0)+i)+ij)).eq.
     :pmap(j,jj))then
      ii=ii+1
      endif
      enddo
      enddo
  319 continue
      elseif(abs(stib(tftyp(0)+i)).eq.3)then
      ii=0
      do 519 j=rhop1,n
      do 515 jj=rdeg(j)+1,rdeg(j)+g(j,j),2
      if(flow(amap(j,jj),0).eq.1)then
      do ij=1,stib(tfnarg(0)+i)
      if(stib(stib(tfo(0)+i)+ij).eq.pmap(j,jj))then
      ii=ii+1
      else
      k=stib(link(0)+stib(stib(tfo(0)+i)+ij))
      if(k.eq.pmap(j,jj))then
      ii=ii+1
      endif
      endif
      enddo
      endif
  515 continue
      do 517 jj=rdeg(j)+g(j,j)+1,degree(j)
      if(flow(amap(j,jj),0).eq.1)then
      do ij=1,stib(tfnarg(0)+i)
      if(stib(stib(tfo(0)+i)+ij).eq.pmap(j,jj))then
      ii=ii+1
      else
      k=stib(link(0)+stib(stib(tfo(0)+i)+ij))
      if(k.eq.pmap(j,jj))then
      ii=ii+1
      endif
      endif
      enddo
      endif
  517 continue
  519 continue
      elseif(abs(stib(tftyp(0)+i)).lt.6)then
      if(abs(stib(tftyp(0)+i)).eq.2)then
      i1=2
      i2=3
      elseif(abs(stib(tftyp(0)+i)).eq.4)then
      i1=2
      i2=2
      elseif(abs(stib(tftyp(0)+i)).eq.5)then
      i1=3
      i2=3
      else
      auxlin(1:srec)='main_3'
      call messag(1,0,0,0)
      endif
      ii=0
      do 419 j=rhop1,n
      do 415 jj=rdeg(j)+1,rdeg(j)+g(j,j),2
      ij=flow(amap(j,jj),0)
      if(ij.eq.i1.or.ij.eq.i2)then
      do ij=1,stib(tfnarg(0)+i)
      if(stib(stib(tfo(0)+i)+ij).eq.pmap(j,jj))then
      ii=ii+1
      else
      k=stib(link(0)+stib(stib(tfo(0)+i)+ij))
      if(k.eq.pmap(j,jj))then
      ii=ii+1
      endif
      endif
      enddo
      endif
  415 continue
      do 417 jj=rdeg(j)+g(j,j)+1,degree(j)
      ij=flow(amap(j,jj),0)
      if(ij.eq.i1.or.ij.eq.i2)then
      do ij=1,stib(tfnarg(0)+i)
      if(stib(stib(tfo(0)+i)+ij).eq.pmap(j,jj))then
      ii=ii+1
      else
      k=stib(link(0)+stib(stib(tfo(0)+i)+ij))
      if(k.eq.pmap(j,jj))then
      ii=ii+1
      endif
      endif
      enddo
      endif
  417 continue
  419 continue
      elseif(abs(stib(tftyp(0)+i)).eq.6)then
      jj=stib(stib(tfo(0)+i)+1)
      ii=0
      i1=stib(vmkvpp(0)+jj)
      i2=stib(vmkvlp(0)+jj)
      do k=rhop1,n
      kk=stib(vfo(k))
      ij=stoz(stcb,stib(i1+kk),stib(i1+kk)-1+stib(i2+kk))
      ii=ii+ij
      enddo
      elseif(abs(stib(tftyp(0)+i)).eq.7)then
      jj=stib(stib(tfo(0)+i)+1)
      ii=0
      i1=stib(pmkvpp(0)+jj)
      i2=stib(pmkvlp(0)+jj)
      do j=rhop1,n
      do k=rdeg(j)+1,rdeg(j)+g(j,j),2
      kk=pmap(j,k)
      ij=stoz(stcb,stib(i1+kk),stib(i1+kk)-1+stib(i2+kk))
      ii=ii+ij
      enddo
      do k=rdeg(j)+g(j,j)+1,degree(j)
      kk=pmap(j,k)
      ij=stoz(stcb,stib(i1+kk),stib(i1+kk)-1+stib(i2+kk))
      ii=ii+ij
      enddo
      enddo
      endif
      if(stib(tftyp(0)+i).gt.0)then
      if(ii.lt.stib(tfa(0)+i).or.ii.gt.stib(tfb(0)+i))then
      goto 58
      endif
      elseif(stib(tftyp(0)+i).lt.0)then
      if(ii.ge.stib(tfa(0)+i).and.ii.le.stib(tfb(0)+i))then
      goto 58
      endif
      endif
  314 continue
      call sdiag(1,0)
      call sdiag(2,0)
      call sdiag(3,0)
      if(jflag(1).eq.0)then
      goto 58
      endif
      if(nloop.eq.0)then
      goto 231
      endif
      do 236 i=rhop1,n
      j=rdeg(i)+1
  420 continue
      if(j.le.degree(i))then
      ii=vmap(i,j)
      aux=g(i,ii)
      k=j+aux
      if(i.ne.ii)then
  430 continue
      if(j.lt.k)then
      kk=1
      aux=pmap(i,j)
  440 continue
      if(j+kk.lt.k)then
      if(aux.eq.pmap(i,j+kk))then
      kk=kk+1
      dsym=dsym*kk
      goto 440
      endif
      endif
      j=j+kk
      goto 430
      endif
      else
  450 continue
      if(j.lt.k)then
      kk=1
      aux=pmap(i,j)
  460 continue
      if(j+kk+kk.lt.k)then
      if(aux.eq.pmap(i,j+kk+kk))then
      kk=kk+1
      dsym=dsym*kk
      goto 460
      endif
      endif
      if(aux.eq.stib(link(0)+aux))then
      do 226 ii=1,kk
      dsym=dsym+dsym
  226 continue
      endif
      j=j+kk+kk
      goto 450
      endif
      endif
      goto 420
      endif
  236 continue
  231 continue
      do i=rhop1,n
      do 250 j=1,degree(i)
      if(vmap(i,j).gt.nleg)then
      k=pmap(i,j)
      if(k.le.stib(link(0)+k))then
      if(k.eq.stib(link(0)+k))then
      if(i.gt.vmap(i,j))then
      goto 250
      endif
      if((i.eq.vmap(i,j)).and.(mod(j-rdeg(i),2).eq.0))then
      goto 250
      endif
      endif
      k=amap(i,j)-nleg
      ex(k)=i
      ey(k)=j
      endif
      endif
  250 continue
      enddo
      do i=rhop1,n
      do j=1,degree(i)
      xli(j)=0
      enddo
      do j=1,degree(i)
      k=stib(stib(vparto(0)+stib(vfo(i)))+j)
      do i3=1,degree(i)
      if(xli(i3).eq.0.and.pmap(i,i3).eq.k)then
      xli(i3)=1
      goto 14
      endif
      enddo
      auxlin(1:srec)='main_4'
      call messag(1,0,0,0)
   14 continue
      ovm(i,j)=i3
      enddo
      enddo
      dis=1
      if(mflag(15).gt.0)then
      goto 495
      endif
      nf=0
      nf1=0
      do i=rhop1,n
      do j=1,degree(i)
      k=ovm(i,j)
      ii=pmap(i,k)
      if(stib(antiq(0)+ii).ne.0)then
      nf=nf+1
      ij=vmap(i,k)
      if(ij.le.nleg)then
      ij=p1(ij)
      if(ij.le.incom)then
      jj=1-2*ij
      else
      jj=2*(incom-ij)
      endif
      elseif(ii.lt.stib(link(0)+ii))then
      jj=2*(amap(i,k)-nleg)-1
      elseif(ii.gt.stib(link(0)+ii))then
      jj=2*(amap(i,k)-nleg)
      elseif(i.lt.ij)then
      jj=2*(amap(i,k)-nleg)-1
      elseif(i.gt.ij)then
      jj=2*(amap(i,k)-nleg)
      elseif(mod(k-rdeg(i),2).ne.0)then
      jj=2*(amap(i,k)-nleg)-1
      else
      jj=2*(amap(i,k)-nleg)
      endif
      if(ij.lt.0)then
      nf1=nf1+1
      endif
      xli(nf)=jj
      endif
      enddo
      enddo
  266 continue
      ii=0
      do i1=1,nf
      if(xli(i1).gt.ii)then
      ii=xli(i1)
      endif
      enddo
      if(ii.gt.0)then
      j1=0
      j2=0
      j3=nf
  293 continue
      if(xli(j3).gt.ii-2)then
      j2=j1
      j1=xli(j3)
      if(j3.ne.nf)then
      xli(j3)=xli(nf)
      dis=-dis
      endif
      nf=nf-1
      endif
      if((j3.gt.1).and.(j2.eq.0))then
      j3=j3-1
      goto 293
      endif
      if(j1.gt.j2)then
      dis=-dis
      endif
      goto 266
      endif
      do i1=1,incom
      do i2=nf,1,-1
      if(xli(i2).eq.1-2*i1)then
      if(i2.ne.nf)then
      xli(i2)=xli(nf)
      dis=-dis
      endif
      nf=nf-1
      goto 292
      endif
      enddo
  292 continue
      enddo
      do i1=rho(1),incom+1,-1
      do i2=nf,1,-1
      if(xli(i2).eq.2*(incom-i1))then
      if(i2.ne.nf)then
      xli(i2)=xli(nf)
      dis=-dis
      endif
      nf=nf-1
      goto 302
      endif
      enddo
  302 continue
      enddo
      if(nf.ne.0)then
      auxlin(1:srec)='main_5'
      call messag(1,0,0,0)
      endif
  495 continue
      call compac
      goto 58
   94 continue
      if(jflag(1).ne.0)then
      call prepos(3)
      endif
      call spp(3)
      stop
      end
      subroutine cbtx(ii,jj,kk)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( sxbuff=2040 )
      character*(srec) auxlin
      common/z22g/auxlin
      character*(sxbuff) stxb
      common/z27g/stxb
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      j1=0
      if(ii.lt.0)then
      j1=1
      elseif(jj.le.0)then
      j1=1
      elseif(ii+jj.gt.srec+ssrec)then
      j1=1
      elseif(kk.lt.0)then
      j1=1
      elseif(kk+jj.gt.scbuff)then
      j1=1
      endif
      if(j1.ne.0)then
      auxlin(1:srec)='cbtx_1'
      call messag(1,0,0,0)
      endif
      stxb(ii+1:ii+jj)=stcb(kk+1:kk+jj)
      ii=ii+jj
      return
      end
      subroutine spp(what)
      implicit integer(a-z)
      save
      parameter ( maxleg=11, maxrho=7, maxi=11, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( sxbuff=2040 )
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1
      common/z3g/nv(0:0),dpntro(0:0),vparto(0:0),vdeg(0:0),nvert
      common/z9g/tpc(0:0)
      common/z11g/nphi,nblok,nprop,npprop
      common/z12g/jflag(1:11),mflag(1:24)
      character*(srec) auxlin
      common/z22g/auxlin
      character*(sxbuff) stxb
      common/z27g/stxb
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z34g/link(0:0),antiq(0:0),blok(0:0)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z44g/xtstrp(0:0),xtstrl(0:0)
      common/z45g/qdatp,qdatl,qvp,qvl,dprefp,dprefl
      common/z47g/ndiagp,ndiagl,hhp,hhl,doffp,doffl,noffp,noffl,wsint
      common/z51g/aplus,aminus,azero,anine,aeq
      integer nps(0:4)
      if(what.eq.0)then
      goto 01
      elseif(what.eq.1)then
      goto 11
      elseif(what.eq.2)then
      goto 21
      elseif(what.eq.3)then
      goto 31
      endif
      auxlin(1:srec)='spp_1'
      call messag(1,0,0,0)
   01 continue
      call hrul(1)
      if(jflag(7).eq.0)then
      print *,''
      endif
      kk=1+(ssrec-qvl)/2
      print *,stxb(1:ssrec)
      print *,stcb(1:kk)//stcb(qvp:qvp-1+qvl)
      print *,stxb(1:ssrec)
      print *,''
      jflag(7)=3
      goto 90
   11 continue
      if(jflag(7).eq.0)then
      print *,''
      endif
      call hrul(1)
      print *,stxb(1:ssrec)
      do i1=1,2
      do i2=0,4
      nps(i2)=0
      enddo
      do i2=1,nphi
      jj=0
      if((i1.eq.1).and.(stib(tpc(0)+i2).eq.5))then
      jj=1
      elseif((i1.eq.2).and.(stib(tpc(0)+i2).ne.5))then
      jj=1
      endif
      if(jj.eq.1)then
      if(i2.lt.stib(link(0)+i2))then
      if(stib(antiq(0)+i2).eq.0)then
      nps(2)=nps(2)+1
      else
      nps(4)=nps(4)+1
      endif
      elseif(i2.eq.stib(link(0)+i2))then
      if(stib(antiq(0)+i2).eq.0)then
      nps(1)=nps(1)+1
      else
      nps(3)=nps(3)+1
      endif
      endif
      endif
      enddo
      nps(0)=nps(1)+nps(2)+nps(3)+nps(4)
      if((i1.eq.2).or.(nps(0).gt.0))then
      ii=0
      call cbtx(ii,2,0)
      call cbtx(ii,stib(xtstrl(0)+i1),stib(xtstrp(0)+i1)-1)
      j1=stcbs(1)
      jj=2
      call vaocb(jj)
      call dkar(nps(0),j1,jj)
      call cbtx(ii,jj,j1)
      call cbtx(ii,stib(xtstrl(0)+14),stib(xtstrp(0)+14)-1)
      if(nps(0).gt.0)then
      do i2=1,4
      if(nps(i2).gt.0)then
      call cbtx(ii,2,0)
      call dkar(nps(i2),j1,jj)
      call cbtx(ii,jj,j1)
      j2=stib(xtstrp(0)+2+i2)
      stcb(j1+1:j1+2)=stcb(j2:j2-1+stib(xtstrl(0)+2+i2))
      call cbtx(ii,2,j1)
      endif
      enddo
      endif
      if(ii.ge.srec)then
      ii=srec-1
      stxb(srec-3:ii)='...'
      endif
      print *,''
      print *,stxb(1:ii)
      endif
      enddo
      ii=0
      call cbtx(ii,2,0)
      call cbtx(ii,stib(xtstrl(0)+7),stib(xtstrp(0)+7)-1)
      call dkar(nvert,j1,jj)
      call cbtx(ii,jj,j1)
      call cbtx(ii,stib(xtstrl(0)+14),stib(xtstrp(0)+14)-1)
      do i1=1,nrho
      if(stib(nv(0)+i1).gt.0)then
      call cbtx(ii,2,0)
      call dkar(i1,j1,jj)
      call cbtx(ii,jj,j1)
      stcb(j1+1:j1+1)='^'
      call cbtx(ii,1,j1)
      call dkar(stib(nv(0)+i1),j1,jj)
      call cbtx(ii,jj,j1)
      endif
      enddo
      if(ii.ge.srec)then
      ii=srec-1
      stxb(srec-3:ii)='...'
      endif
      print *,''
      print *,stxb(1:ii)
      print *,''
      call hrul(1)
      print *,stxb(1:ssrec)
      print *,''
      jflag(7)=3
      goto 90
   21 continue
      gg=0
      do i1=3,nrho
      gg=gg+(i1-2)*rho(i1)
      enddo
      call vaocb(1)
      j1=stcbs(1)
      j2=j1+1
      ii=0
      do i1=1,nrho
      if(stib(nv(0)+i1).gt.0)then
      kk=wztos(gg/(i1-2))+2
      call dkar(i1,j1,jj)
      if(rho(i1).gt.0)then
      call cbtx(ii,jj,j1)
      stcb(j2:j2)='^'
      call cbtx(ii,1,j1)
      call dkar(rho(i1),j1,jj)
      call cbtx(ii,jj,j1)
      kk=kk-jj
      else
      call cbtx(ii,jj,0)
      stcb(j2:j2)=char(aminus)
      call cbtx(ii,1,j1)
      endif
      call cbtx(ii,kk,0)
      endif
      enddo
      ii=ii-2
      j2=ii
      jj=5
      call cbtx(ii,jj,0)
      call cbtx(ii,stib(xtstrl(0)+15),stib(xtstrp(0)+15)-1)
      call cbtx(ii,jj,0)
      call cbtx(ii,hhl,hhp)
      if(jflag(3).ne.0)then
      call cbtx(ii,stib(xtstrl(0)+16),stib(xtstrp(0)+16)-1)
      endif
      i2=stib(xtstrl(0)+8)
      j1=stib(xtstrl(0)+15)
      if(j2.lt.i2)then
      j4=(i2-j2)/2
      else
      j4=0
      endif
      if(jflag(8).eq.0)then
      if(j2.lt.i2)then
      j3=0
      else
      j3=(j2-i2)/2
      endif
      jj=j2-i2+2*jj+j1+j4-j3-3
      i1=stib(xtstrp(0)+8)
      i3=stib(xtstrp(0)+9+mflag(11))
      i4=i3-1+stib(xtstrl(0)+9+mflag(11))
      print *,stcb(1:j3+4)//stcb(i1:i1-1+i2)//stcb(1:jj)//stcb(i3:i4)
      print *,''
      jflag(8)=3
      endif
      print *,stcb(1:j4+4)//stxb(1:ii)
      jflag(7)=0
      goto 90
   31 continue
      ii=0
      call cbtx(ii,stib(xtstrl(0)+11),stib(xtstrp(0)+11)-1)
      call cbtx(ii,ndiagl,ndiagp)
      jj=0
      if(ndiagl.eq.1)then
      j1=ndiagp+1
      if(ichar(stcb(j1:j1)).eq.azero+1)then
      jj=1
      endif
      endif
      j1=12+mflag(11)
      call cbtx(ii,stib(xtstrl(0)+j1)-jj,stib(xtstrp(0)+j1)-1)
      if(jflag(2).ne.0)then
      call cbtx(ii,stib(xtstrl(0)+16),stib(xtstrp(0)+16)-1)
      endif
      if(jflag(7).eq.0)then
      print *,''
      endif
      print *,''
      print *,stxb(1:ii)
      print *,''
      jflag(7)=3
   90 continue
      return
      end
      subroutine rflag
      implicit integer(a-z)
      save
      parameter ( maxleg=11, maxrho=7, maxi=11, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      common/z1in/leg(1:maxleg),nleg,nloop,incom
      common/z9g/tpc(0:0)
      common/z12g/jflag(1:11),mflag(1:24)
      common/z14g/zcho(0:maxli),zbri(0:maxli),zpro(0:maxli),
     :rbri(0:maxli),sbri(0:maxli)
      character*(srec) auxlin
      common/z22g/auxlin
      common/z30g/stib(1:sibuff)
      if(mflag(11).ne.0)then
      if(mflag(14).eq.0)then
      auxlin(1:srec)='option "topol" does not apply here'
      call messag(1,0,0,0)
      endif
      ii=0
      if(mflag(22).ne.0)then
      ii=1
      elseif(mflag(23).ne.0)then
      ii=1
      endif
      if(ii.ne.0)then
      auxlin(1:srec)='elink, plink incompatible with option "topol"'
      call messag(1,0,0,0)
      endif
      endif
      if(mflag(12).ne.0)then
      if((mflag(15).gt.0).or.(mflag(16).eq.0))then
      auxlin(1:srec)='option "floop" does not apply here'
      call messag(1,0,0,0)
      endif
      endif
      if(mflag(1).gt.0)then
      if(nleg.ne.1)then
      if(mflag(2).lt.0)then
      jflag(2)=1
      elseif(mflag(3).lt.0)then
      jflag(2)=1
      elseif(mflag(4).lt.0)then
      jflag(2)=1
      endif
      endif
      endif
      if(mflag(6).gt.0)then
      if(mflag(4).lt.0)then
      jflag(2)=1
      endif
      endif
      if(nloop.eq.0)then
      do i1=1,maxli
      zcho(i1)=0
      enddo
      endif
      if(mflag(1).gt.0)then
      do i1=1,maxli
      zbri(i1)=0
      enddo
      elseif(mflag(1).lt.0)then
      zbri(0)=0
      endif
      if(mflag(2).gt.0)then
      do i1=1,maxli
      sbri(i1)=0
      enddo
      elseif(mflag(2).lt.0)then
      sbri(0)=0
      zcho(0)=0
      endif
      if(mflag(4).lt.0)then
      zbri(0)=0
      zcho(0)=0
      endif
      if(mflag(5).lt.0)then
      zcho(0)=0
      endif
      if(mflag(6).lt.0)then
      zcho(0)=0
      endif
      if(mflag(8).gt.0)then
      do i1=1,maxli
      sbri(i1)=0
      enddo
      elseif(mflag(8).lt.0)then
      zcho(0)=0
      endif
      if(mflag(9).lt.0)then
      zcho(0)=0
      endif
      if(mflag(3).ne.0)then
      jflag(10)=1
      endif
      if(mflag(5).ne.0)then
      jflag(10)=1
      endif
      if(mflag(23).ne.0)then
      jflag(10)=1
      endif
      if(nloop.eq.0)then
      if(mflag(2).ge.0)then
      mflag(2)=0
      else
      jflag(2)=1
      endif
      if(mflag(4).ge.0)then
      mflag(4)=0
      else
      jflag(2)=1
      endif
      if(mflag(5).ge.0)then
      mflag(5)=0
      else
      jflag(2)=1
      endif
      if(mflag(6).ge.0)then
      mflag(6)=0
      else
      jflag(2)=1
      endif
      if(mflag(8).ge.0)then
      mflag(8)=0
      else
      jflag(2)=1
      endif
      if(mflag(9).ge.0)then
      mflag(9)=0
      else
      jflag(2)=1
      endif
      if(mflag(12).ge.0)then
      mflag(12)=0
      else
      jflag(2)=1
      endif
      endif
      if(nloop.eq.0)then
      mflag(19)=0
      elseif(mflag(1).gt.0)then
      mflag(19)=0
      elseif(mflag(2).gt.0)then
      mflag(19)=0
      elseif(mflag(8).gt.0)then
      mflag(19)=0
      endif
      if(nleg.eq.1)then
      if(stib(tpc(0)+leg(1)).eq.1)then
      jflag(2)=1
      endif
      endif
      return
      end
      subroutine prepos(what)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( eoa=-127, nap=-2047 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      common/z5in/dunit,munit,ounit,sunit,funit
      common/z9in/lofile,ofilea,ofileb
      common/z24g/iogp(1:4)
      character*(srec) auxlin
      common/z22g/auxlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      common/z37g/drecp(0:0),drecl(0:0),drecii(0:0),irecc(0:0),
     :frecc(0:0),ndrec,ncom
      common/z38g/gmkp(0:0),gmkl(0:0),gmkd(0:0),gmko(0:0),
     :gmkvp(0:0),gmkvl(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      common/z45g/qdatp,qdatl,qvp,qvl,dprefp,dprefl
      common/z47g/ndiagp,ndiagl,hhp,hhl,doffp,doffl,noffp,noffl,wsint
      common/z52g/acf0(0:127),acf1(0:127),lfeed,aspace,dquote,squote
      if((what.ne.1).and.(what.ne.3))then
      auxlin(1:srec)='prepos_1'
      call messag(1,0,0,0)
      endif
      if(what.eq.1)then
      ounit=5
      call iopen(ofilea,lofile,ounit,1)
      endif
      ig=iogp(what)
      klin=0
      lupt=0
   20 continue
      if((ig.lt.iogp(what)).or.(ig.ge.iogp(what+1)))then
      goto 80
      endif
      if(stib(ig).gt.0)then
      klin=klin+1
      call vaocb(klin)
      i=stcbs(1)+klin
      stcb(i:i)=char(stib(ig))
      ig=ig+1
      goto 20
      endif
      if(stib(ig).eq.-1)then
      if(stib(ig+1).eq.11)then
      if(stib(ig+4).eq.0)then
      stib(ig+4)=ncom
      stib(ig+5)=0
      endif
      stib(ig+5)=stib(ig+5)+1
      if(stib(ig+4).lt.stib(ig+5))then
      stib(ig+4)=0
      stib(ig+5)=0
      ig=ig+1
      else
      lupt=11
      lupi=stib(ig+5)
      endif
      elseif(stib(ig+1).eq.12)then
      if(stib(ig+4).eq.0)then
      stib(ig+4)=stib(frecc(0)+lupi)
      stib(ig+5)=stib(irecc(0)+lupi)-1
      endif
      stib(ig+5)=stib(ig+5)+1
      if(stib(ig+4).lt.stib(ig+5))then
      stib(ig+4)=0
      stib(ig+5)=0
      ig=ig+1
      else
      lupt=12
      lupj=stib(ig+5)
      endif
      elseif(stib(ig+1).eq.19)then
      if(lupt.eq.11)then
      lupt=0
      else
      lupt=11
      endif
      else
      goto 80
      endif
      elseif(stib(ig).eq.-2)then
      if(stib(ig+1).eq.21)then
      if(klin.le.0)then
      goto 80
      endif
      klin=klin-1
      else
      goto 80
      endif
      elseif(stib(ig).eq.-3)then
      if(stib(ig+1).eq.71)then
      i=stcbs(1)+klin
      klin=klin+ndiagl
      call vaocb(klin)
      stcb(i+1:i+ndiagl)=stcb(ndiagp+1:ndiagp+ndiagl)
      elseif(stib(ig+1).eq.81)then
      klin=klin+qvl
      call vaocb(klin)
      i=stcbs(1)+klin
      stcb(i-qvl+1:i)=stcb(qvp:qvp-1+qvl)
      else
      goto 80
      endif
      elseif(stib(ig).eq.-4)then
      if(stib(ig+1).eq.82)then
      if(lupt.eq.11)then
      do 70 i1=stib(irecc(0)+lupi),stib(frecc(0)+lupi)
      j=stib(drecp(0)+i1)
      jj=stib(drecl(0)+i1)
      klin=klin+jj
      call vaocb(klin)
      i=stcbs(1)+klin
      stcb(i-jj+1:i)=stcb(j:j-1+jj)
   70 continue
      elseif(lupt.eq.12)then
      jj=stib(drecl(0)+lupj)
      klin=klin+jj
      call vaocb(klin)
      j=stib(drecp(0)+lupj)
      i=stcbs(1)+klin
      stcb(i-jj+1:i)=stcb(j:j-1+jj)
      else
      goto 80
      endif
      else
      goto 80
      endif
      elseif(stib(ig).eq.-6)then
      ii=stib(ig+1)
      if((ii.le.0).or.(ii.gt.nudk))then
      goto 80
      endif
      if(stib(udkt(0)+ii).eq.1)then
      ij=stib(udki(0)+ii)
      if((ij.le.0).or.(ij.gt.ngmk))then
      goto 80
      elseif(stib(gmkd(0)+ij).ne.1)then
      goto 80
      endif
      ij=stib(gmko(0)+ij)+1
      j=stib(gmkvp(0)+ij)
      jj=stib(gmkvl(0)+ij)
      if(jj.gt.0)then
      klin=klin+jj
      call vaocb(klin)
      i=stcbs(1)+klin
      stcb(i-jj+1:i)=stcb(j:j-1+jj)
      endif
      if(jj.lt.0)then
      goto 80
      endif
      else
      goto 80
      endif
      elseif(stib(ig).eq.eoa)then
      goto 30
      else
      goto 80
      endif
      ig=ig+stib(ig+2)
      goto 20
   30 continue
      if(klin.gt.0)then
      i=stcbs(1)+klin
      if(ichar(stcb(i:i)).ne.lfeed)then
      goto 80
      endif
      endif
      j1=stcbs(1)+1
      do 29 i1=stcbs(1)+1,stcbs(1)+klin
      if(ichar(stcb(i1:i1)).eq.lfeed)then
      ios=1
      if(j1.lt.i1)then
      write(unit=ounit,fmt=404,iostat=ios)stcb(j1:i1-1)
      else
      write(unit=ounit,fmt=404,iostat=ios)
      endif
      if(ios.ne.0)then
      auxlin(1:srec)='run-time error while writing to'
      call messag(1,0,0,5)
      endif
      j1=i1+1
      endif
   29 continue
  404 format(a)
      if(what.eq.3)then
      ios=1
      close(unit=ounit,status='keep',iostat=ios)
      if(ios.ne.0)then
      auxlin(1:srec)='system/filesystem error'
      call messag(1,0,0,0)
      endif
      ounit=0
      endif
      goto 90
   80 continue
      if(what.eq.1)then
      auxlin(1:srec)='run-time error while processing <prologue>'
      else
      auxlin(1:srec)='run-time error while processing <epilogue>'
      endif
      call messag(1,0,0,0)
   90 continue
      return
      end
      subroutine umpi(xx,situ)
      implicit integer(a-z)
      save
      parameter ( maxleg=11, maxrho=7, maxi=11, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( srec=81, ssrec=62 )
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1
      common/z4g/n,nli
      common/z17g/xtail(1:maxn),xhead(1:maxn),ntadp
      character*(srec) auxlin
      common/z22g/auxlin
      integer aa(1:maxn),bb(1:maxn)
      ntadp=0
      situ=-1
      do 66 ii=rhop1,n
      do 56 jj=ii+1,n
      if(g(ii,jj).eq.1)then
      g(ii,jj)=0
      k=1
      bb(1)=1
      aa(1)=1
      do 06 i=2,n
      aa(i)=0
   06 continue
      do 46 j=1,n
      if(k.eq.n)then
      goto 20
      endif
      if(j.gt.k)then
      g(ii,jj)=1
      if(xx.eq.1)then
      goto 90
      endif
      aux=0
      do 16 kk=1,rho(1)
      aux=aux+aa(kk)
   16 continue
      if(xx.eq.2)then
      if(aux.eq.0.or.aux.eq.rho(1))then
      goto 90
      endif
      elseif(xx.eq.3)then
      if(aux.eq.0.or.aux.eq.rho(1))then
      ntadp=ntadp+1
      xtail(ntadp)=ii
      xhead(ntadp)=jj
      endif
      elseif(xx.eq.4)then
      if(aux.eq.1.or.aux.eq.rho(1)-1)then
      goto 90
      endif
      endif
      goto 20
      endif
      aux=bb(j)
      do 26 i=1,aux-1
      if(aa(i).eq.0)then
      if(g(i,aux).gt.0)then
      k=k+1
      bb(k)=i
      aa(i)=1
      endif
      endif
   26 continue
      do 36 i=aux+1,n
      if(aa(i).eq.0)then
      if(g(aux,i).gt.0)then
      k=k+1
      bb(k)=i
      aa(i)=1
      endif
      endif
   36 continue
   46 continue
      auxlin(1:srec)='umpi_1'
      call messag(1,0,0,0)
   20 continue
      g(ii,jj)=1
      endif
   56 continue
   66 continue
      situ=1
   90 continue
      return
      end
      subroutine umvi(xx,situ)
      implicit integer(a-z)
      save
      parameter ( maxleg=11, maxrho=7, maxi=11, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( srec=81, ssrec=62 )
      common/z1in/leg(1:maxleg),nleg,nloop,incom
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1
      common/z4g/n,nli
      character*(srec) auxlin
      common/z22g/auxlin
      integer aa(1:maxn),bb(1:maxn),cc(1:maxn)
      if(xx.eq.1)then
      if(nloop.eq.0)then
      goto 70
      endif
      ii=0
      if(rho(1).ne.1)then
      ii=1
      elseif(nloop.gt.1)then
      ii=1
      endif
      if(ii.ne.0)then
      do i1=rhop1,n
      if(g(i1,i1).ne.0)then
      goto 80
      endif
      enddo
      endif
      elseif(xx.eq.2)then
      if(nloop.eq.0)then
      goto 70
      elseif(rho(1).eq.0)then
      goto 70
      endif
      goto 200
      else
      auxlin(1:srec)='umvi_0'
      call messag(1,0,0,0)
      endif
      do 77 i1=rhop1,n
      kk=rho(1)
      do 02 i2=1,rho(1)
      aa(i2)=1
      bb(i2)=i2
   02 continue
      do 04 i2=rhop1,n
      aa(i2)=0
   04 continue
      do 46 i2=1,n
      if(kk.eq.n)then
      goto 77
      endif
      if(i2.gt.kk)then
      goto 80
      endif
      j1=bb(i2)
      if(j1.ne.i1)then
      do i3=rhop1,j1-1
      if(aa(i3).eq.0)then
      if(g(i3,j1).gt.0)then
      kk=kk+1
      bb(kk)=i3
      aa(i3)=1
      endif
      endif
      enddo
      do i3=j1+1,n
      if(aa(i3).eq.0)then
      if(g(j1,i3).gt.0)then
      kk=kk+1
      bb(kk)=i3
      aa(i3)=1
      endif
      endif
      enddo
      endif
   46 continue
      auxlin(1:srec)='umvi_1'
      call messag(1,0,0,0)
   77 continue
      goto 70
  200 continue
      do 277 i1=rhop1,n
      do i2=1,n
      aa(i2)=0
      enddo
      jj=0
      kk=1
      aa(i1)=-1
      bb(1)=i1
      do 146 i2=1,n
      if(kk.eq.n)then
      goto 150
      elseif(aa(i2).ne.0)then
      goto 146
      endif
      jj=jj+1
      aa(i2)=jj
      kk=kk+1
      bb(kk)=i2
      ii=kk
  106 continue
      j1=bb(ii)
      do i3=i2+1,j1-1
      if(aa(i3).eq.0)then
      if(g(i3,j1).gt.0)then
      kk=kk+1
      bb(kk)=i3
      aa(i3)=jj
      endif
      endif
      enddo
      do i3=j1+1,n
      if(aa(i3).eq.0)then
      if(g(j1,i3).gt.0)then
      kk=kk+1
      bb(kk)=i3
      aa(i3)=jj
      endif
      endif
      enddo
      if(ii.lt.kk)then
      if(kk.lt.n)then
      ii=ii+1
      goto 106
      endif
      endif
  146 continue
  150 continue
      if(kk.ne.n)then
      auxlin(1:srec)='umvi_2'
      call messag(1,0,0,0)
      endif
      do i2=1,jj
      cc(i2)=0
      enddo
      do i2=1,rho(1)
      cc(aa(i2))=cc(aa(i2))+1
      enddo
      do i2=1,jj
      if(cc(i2).eq.1)then
      goto 160
      endif
      enddo
      goto 277
  160 continue
      j0=0
      do i2=1,jj
      bb(i2)=cc(i2)
      if(cc(i2).eq.0)then
      j0=j0+1
      endif
      enddo
      do i2=rhop1,n
      if(i2.ne.i1)then
      bb(aa(i2))=bb(aa(i2))+1
      endif
      enddo
      do 170 i2=1,jj
      if(cc(i2).eq.1)then
      j1=bb(i2)
      j2=g(i1,i1)
      j3=j0
      if(j1.eq.1)then
      if(j2.gt.1)then
      j2=j2-2
      elseif(j0.gt.0)then
      i4=0
      do i3=1,jj
      if(cc(i3).eq.0)then
      if(i4.eq.0)then
      i4=i3
      elseif(bb(i3).lt.bb(i4))then
      i4=i3
      endif
      endif
      enddo
      j3=j3-1
      j1=j1+bb(i4)
      else
      goto 170
      endif
      endif
      if(0.lt.j2+j3)then
      goto 80
      elseif(j1+2.lt.n)then
      goto 80
      endif
      endif
  170 continue
  277 continue
   70 continue
      situ=1
      goto 90
   80 continue
      situ=-1
   90 continue
      return
      end
      subroutine dkar(im,ic,iw)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      character*(srec) auxlin
      common/z22g/auxlin
      character*(scbuff) stcb
      common/z31g/stcb
      common/z48g/iref,wiref
      common/z51g/aplus,aminus,azero,anine,aeq
      character*(1) c1
      jc=ic
      if(im.lt.0)then
      if(im.lt.-iref)then
      goto 30
      endif
      xm=-im
      if(jc.ge.scbuff)then
      iw=wiref
      call vaocb(iw)
      endif
      jc=jc+1
      stcb(jc:jc)=char(aminus)
      else
      if(im.gt.iref)then
      goto 30
      endif
      xm=im
      endif
      ii=jc+1
   10 continue
      if(jc.ge.scbuff)then
      iw=(jc-ic)+wiref
      call vaocb(iw)
      endif
      jc=jc+1
      if(xm.lt.10)then
      stcb(jc:jc)=char(xm+azero)
      else
      ym=xm/10
      stcb(jc:jc)=char((xm-10*ym)+azero)
      xm=ym
      goto 10
      endif
      jj=jc
   20 continue
      if(ii.lt.jj)then
      c1(1:1)=stcb(ii:ii)
      stcb(ii:ii)=stcb(jj:jj)
      stcb(jj:jj)=c1(1:1)
      ii=ii+1
      jj=jj-1
      goto 20
      endif
      goto 90
   30 continue
      auxlin(1:srec)='dkar_1'
      call messag(1,0,0,0)
   90 continue
      iw=jc-ic
      return
      end
      subroutine sdiag(ii,jj)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      character*(srec) auxlin
      common/z22g/auxlin
      character*(scbuff) stcb
      common/z31g/stcb
      common/z47g/ndiagp,ndiagl,hhp,hhl,doffp,doffl,noffp,noffl,wsint
      common/z51g/aplus,aminus,azero,anine,aeq
      if(ii.eq.1)then
      xp=hhp
      xl=hhl
      elseif(ii.eq.2)then
      xp=noffp
      xl=noffl
      elseif(ii.eq.3)then
      xp=ndiagp
      xl=ndiagl
      else
      auxlin(1:srec)='sdiag_1'
      call messag(1,0,0,0)
      endif
      if(jj.eq.0)then
      do i1=xp+xl,xp+1,-1
      j1=ichar(stcb(i1:i1))
      if(j1.lt.anine)then
      stcb(i1:i1)=char(j1+1)
      goto 90
      else
      stcb(i1:i1)=char(azero)
      endif
      enddo
      if(xl.ge.wsint)then
      auxlin(1:srec)='sdiag_2'
      call messag(1,0,0,0)
      endif
      j1=xp+1
      stcb(j1:j1)=char(azero+1)
      j1=j1+xl
      stcb(j1:j1)=char(azero)
      xl=xl+1
      elseif(jj.eq.-1)then
      j1=xp+1
      stcb(j1:j1)=char(azero)
      xl=1
      else
      auxlin(1:srec)='sdiag_3'
      call messag(1,0,0,0)
      endif
      if(ii.eq.1)then
      hhl=xl
      elseif(ii.eq.2)then
      noffl=xl
      elseif(ii.eq.3)then
      ndiagl=xl
      endif
   90 continue
      return
      end
      integer function stoz(s,j1,j2)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      character*(*) s
      character*(srec) auxlin
      common/z22g/auxlin
      common/z48g/iref,wiref
      common/z51g/aplus,aminus,azero,anine,aeq
      common/z52g/acf0(0:127),acf1(0:127),lfeed,aspace,dquote,squote
      stoz=0
      if((j1.le.0).or.(j1.gt.j2))then
      goto 30
      endif
      is=ichar(s(j1:j1))
      if((is.eq.aplus).or.(is.eq.aminus))then
      i1=j1+1
      else
      i1=j1
      endif
   10 continue
      jj=ichar(s(i1:i1))
      if(acf1(jj).ne.1)then
      goto 30
      endif
      stoz=10*stoz+(jj-azero)
      if(abs(stoz).gt.iref)then
      call uput(4)
      endif
      if(i1.lt.j2)then
      i1=i1+1
      goto 10
      endif
      if(is.eq.aminus)then
      stoz=-stoz
      endif
      goto 90
   30 continue
      auxlin(1:srec)='stoz_1'
      call messag(1,0,0,0)
   90 continue
      return
      end
      integer function nullz(s,ia,ib)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      character*(*) s
      character*(srec) auxlin
      common/z22g/auxlin
      common/z51g/aplus,aminus,azero,anine,aeq
      common/z52g/acf0(0:127),acf1(0:127),lfeed,aspace,dquote,squote
      if(ia.le.0)then
      goto 30
      endif
      nullz=0
      jj=ichar(s(ia:ia))
      if((jj.eq.aplus).or.(jj.eq.aminus))then
      ii=1
      else
      ii=0
      endif
      if(ia+ii.gt.ib)then
      goto 30
      endif
      do i1=ia+ii,ib
      jj=ichar(s(i1:i1))
      if(acf1(jj).ne.1)then
      goto 30
      elseif(jj.ne.azero)then
      goto 90
      endif
      enddo
      nullz=1
      goto 90
   30 continue
      auxlin(1:srec)='nullz_1'
      call messag(1,0,0,0)
   90 continue
      return
      end
      integer function wztos(nn)
      implicit integer(a-z)
      save
      if(nn.lt.0)then
      mm=-nn
      wztos=2
      else
      mm=nn
      wztos=1
      endif
   10 continue
      if(mm.gt.9)then
      if(mm.lt.100)then
      wztos=wztos+1
      else
      mm=mm/100
      wztos=wztos+2
      goto 10
      endif
      endif
      return
      end
      subroutine spak(s,ind,m,uc,nos)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      character*(srec) auxlin
      common/z22g/auxlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z51g/aplus,aminus,azero,anine,aeq
      common/z52g/acf0(0:127),acf1(0:127),lfeed,aspace,dquote,squote
      character*(*) s
      if((cmal(0).ne.1).or.(imal(0).ne.1))then
      goto 08
      elseif((m.ne.0).and.(m.ne.1))then
      goto 08
      elseif((nos.ne.0).and.(nos.ne.1))then
      goto 08
      elseif((uc.ne.0).and.(uc.ne.1))then
      goto 08
      endif
      inos=nos
      iuc=uc
      i1=0
      sl=1
      ii=2*srec+2
   03 continue
      if(s(sl:sl).ne.';')then
      if(sl.lt.ii)then
      if((i1.eq.0).and.(s(sl:sl).eq.','))then
      i1=sl
      endif
      sl=sl+1
      goto 03
      else
      goto 08
      endif
      endif
      if(i1.eq.0)then
      i1=sl
      endif
      if((i1.lt.3).or.(sl.lt.3))then
      goto 08
      elseif(mod(i1,2).eq.0)then
      goto 08
      endif
      sp=stcbs(1)+1
      i1=(i1-1)/2
      call aocb(i1+1)
      ij=0
      j=1
      k=sp
      do i2=1,i1
      ii=ichar(s(j:j))
      jj=ichar(s(j+1:j+1))
      if((ii.lt.azero).or.(ii.gt.anine))then
      goto 08
      elseif((jj.lt.azero).or.(jj.gt.anine))then
      goto 08
      endif
      kk=10*ii+jj-496
      if((iuc.eq.1).and.(kk.gt.64).and.(kk.lt.91))then
      goto 08
      endif
      if((kk.gt.96).and.(kk.lt.123))then
      ij=1
      endif
      stcb(k:k)=char(kk)
      j=j+2
      k=k+1
      enddo
      stcb(k:k)=char(lfeed)
      if((uic.eq.1).and.(ij.eq.0))then
      goto 08
      endif
      if(inos.eq.1)then
      ichs=115
      if(ichar(stcb(k-1:k-1)).ne.ichs)then
      inos=0
      endif
      endif
      j0=stibs(1)
      if(m.eq.0)then
      call aoib(2)
      elseif(m.eq.1)then
      call vaoib(2)
      endif
      stib(j0+1)=sp
      stib(j0+2)=i1
      if(m.eq.1)then
      if((nos.ne.0).or.(uc.ne.0))then
      goto 08
      endif
      goto 90
      endif
      j2=i1+i1+1
      j1=j2+1
   11 continue
      if(j2.lt.sl)then
      j2=j2+1
      if((s(j2:j2).eq.',').or.(s(j2:j2).eq.';'))then
      j2=j2-1
      if(j2.lt.j1)then
      goto 08
      endif
      i3=stoz(s,j1,j2)
      call aoib(1)
      stib(stibs(1))=i3
      j2=j2+1
      j1=j2+1
      endif
      goto 11
      endif
      ind=ind+1
      if(iuc.eq.1)then
      ii=stcbs(1)
      call aocb(i1+1)
      i3=ii+1
      i4=ii-i1
      sp=i3
      do i2=1,i1
      j=ichar(stcb(i4:i4))
      if((j.gt.96).and.(j.lt.123))then
      j=j-32
      endif
      stcb(i3:i3)=char(j)
      i3=i3+1
      i4=i4+1
      enddo
      stcb(i3:i3)=stcb(i4:i4)
      j1=stibs(1)
      jj=j1-j0
      if(jj.lt.2)then
      goto 08
      endif
      call aoib(jj)
      ii=j1+1
      stib(ii)=stib(ii-jj)+i1+1
      do i2=2,jj
      ii=ii+1
      stib(ii)=stib(ii-jj)
      enddo
      ind=ind+1
      endif
      if(inos.eq.1)then
      j1=stibs(1)
      jj=j1-j0
      if(jj.lt.4)then
      goto 08
      endif
      call aoib(jj)
      ii=j1
      do i2=1,jj
      ii=ii+1
      stib(ii)=stib(ii-jj)
      enddo
      ii=j1+2
      stib(ii)=stib(ii)-1
      ind=ind+1
      if(uc.eq.1)then
      ii=ii+jj/2
      stib(ii)=stib(ii)-1
      ind=ind+1
      endif
      endif
      goto 90
   08 continue
      auxlin(1:srec)='spak_1'
      call messag(1,0,0,0)
   90 continue
      return
      end
      integer function stds(s,ia,ib,iw,ilf)
      implicit integer(a-z)
      save
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z52g/acf0(0:127),acf1(0:127),lfeed,aspace,dquote,squote
      character*(*) s
      stds=-1
      if((ia.le.0).or.(ia.ge.ib))then
      goto 90
      elseif((iw.lt.0).or.(iw.gt.2))then
      goto 90
      elseif((ilf.ne.0).and.(ilf.ne.1))then
      goto 90
      endif
      sc=0
      sl=0
      sf=0
      itl=1
      ii=0
      ll=0
      mm=0
      do i1=ia,ib
      jj=ichar(s(i1:i1))
      if(jj.eq.squote)then
      ii=ii+1
      mm=1-mm
      ll=0
      if(ii.gt.1)then
      kk=mm
      else
      kk=0
      endif
      elseif(jj.eq.lfeed)then
      if((mm.ne.0).or.(ll.eq.1))then
      goto 90
      endif
      sf=sf+1
      ll=1
      if(ii.gt.0)then
      ii=0
      sc=sc+1
      endif
      kk=0
      elseif(jj.eq.aspace)then
      if(ii.gt.0)then
      ii=0
      if(mm.eq.0)then
      sc=sc+1
      endif
      endif
      kk=mm
      else
      if((mm.ne.1).or.(acf0(jj).lt.0))then
      goto 90
      endif
      ii=0
      kk=1
      endif
      if(kk.ne.0)then
      if(iw.eq.2)then
      if(jj.ne.aspace)then
      itl=0
      endif
      endif
      if((iw.eq.1).or.(iw.eq.2.and.itl.eq.0))then
      sl=sl+1
      call vaocb(sl)
      j1=stcbs(1)+sl
      stcb(j1:j1)=s(i1:i1)
      endif
      endif
      enddo
      if(mm.ne.0)then
      goto 90
      endif
      if(ichar(s(ib:ib)).eq.squote)then
      sc=sc+1
      endif
      if((iw.eq.2).and.(sl.gt.1))then
      do i1=stcbs(1)+sl,stcbs(1)+1,-1
      if(ichar(stcb(i1:i1)).eq.aspace)then
      sl=sl-1
      else
      goto 20
      endif
      enddo
      endif
   20 continue
      stds=sl
   90 continue
      return
      end
      subroutine mstr0(s,ia,ib,pp,pl,ind)
      implicit integer(a-z)
      save
      parameter ( eoa=-127, nap=-2047 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      character*(*) s
      ind=0
      if(pp.eq.nap)then
      goto 90
      endif
      ls=ib-ia+1
   10 continue
      ind=ind+1
      ii=stib(pl+ind)
      if(ii.ne.eoa)then
      if(ii.eq.ls)then
      jj=stib(pp+ind)
      if(stcb(jj:jj-1+ls).eq.s(ia:ib))then
      goto 90
      endif
      endif
      goto 10
      endif
      ind=0
   90 continue
      return
      end
      subroutine mstr1(s,ia,ib,pp,kp,kl,ind)
      implicit integer(a-z)
      save
      parameter ( eoa=-127, nap=-2047 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      character*(*) s
      ind=0
      if(pp.eq.nap)then
      goto 90
      endif
      ls=ib-ia+1
      p1=pp
   10 continue
      p1=stib(p1)
      if(p1.eq.eoa)then
      ind=0
      goto 90
      endif
      ind=ind+1
      if(stib(p1+kl).eq.ls)then
      jj=stib(p1+kp)
      if(stcb(jj:jj-1+ls).eq.s(ia:ib))then
      goto 90
      endif
      endif
      goto 10
   90 continue
      return
      end
      integer function stdw(s,ia,ib)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      character*(srec) auxlin
      common/z22g/auxlin
      common/z52g/acf0(0:127),acf1(0:127),lfeed,aspace,dquote,squote
      character*(*) s
      if(ia.gt.ib.or.ia.le.0)then
      auxlin(1:srec)='stdw_1'
      call messag(1,0,0,0)
      endif
      if(acf1(ichar(s(ia:ia))).lt.2)then
      stdw=0
      else
      stdw=1
      endif
      if(stdw.eq.1)then
      do i1=ia+1,ib
      if(acf1(ichar(s(i1:i1))).lt.0)then
      stdw=0
      goto 90
      endif
      enddo
      endif
   90 continue
      return
      end
      integer function stdz(s,ia,ib)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      character*(srec) auxlin
      common/z22g/auxlin
      common/z51g/aplus,aminus,azero,anine,aeq
      common/z52g/acf0(0:127),acf1(0:127),lfeed,aspace,dquote,squote
      character*(*) s
      if((ia.le.0).or.(ia.gt.ib))then
      auxlin(1:srec)='stdz_1'
      call messag(1,0,0,0)
      endif
      stdz=0
      jj=ichar(s(ia:ia))
      if((jj.eq.aplus).or.(jj.eq.aminus))then
      if(ia.lt.ib)then
      stdz=1
      endif
      elseif(acf1(jj).eq.1)then
      stdz=1
      endif
      if(stdz.eq.1)then
      do i1=ia+1,ib
      if(acf1(ichar(s(i1:i1))).ne.1)then
      stdz=0
      goto 90
      endif
      enddo
      endif
   90 continue
      return
      end
      integer function stdq(s,ia,ib)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      character*(srec) auxlin
      common/z22g/auxlin
      common/z51g/aplus,aminus,azero,anine,aeq
      common/z52g/acf0(0:127),acf1(0:127),lfeed,aspace,dquote,squote
      character*(*) s
      if((ia.le.0).or.(ia.gt.ib))then
      auxlin(1:srec)='stdq_1'
      call messag(1,0,0,0)
      endif
      xsla=0
      xnum=0
      xden=0
      jj=ichar(s(ia:ia))
      if(jj.eq.aplus)then
      stdq=1
      elseif(jj.eq.aminus)then
      stdq=1
      elseif(acf1(jj).eq.1)then
      stdq=1
      xnum=1
      else
      stdq=0
      endif
      if(stdq.eq.1)then
      do i1=ia+1,ib
      jj=ichar(s(i1:i1))
      if(acf1(jj).eq.1)then
      if(xsla.eq.0)then
      xnum=1
      elseif(jj.ne.azero)then
      xden=1
      endif
      elseif(s(i1:i1).eq.'/')then
      if(xsla.ne.0)then
      stdq=0
      goto 90
      endif
      xsla=1
      else
      stdq=0
      goto 90
      endif
      enddo
      if(xnum.ne.1)then
      stdq=0
      elseif(xsla.ne.0)then
      if(xden.ne.1)then
      stdq=0
      endif
      endif
      endif
   90 continue
      return
      end
      subroutine compac
      implicit integer(a-z)
      save
      parameter ( maxleg=11, maxrho=7, maxi=11, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( srec=81, ssrec=62 )
      parameter ( eoa=-127, nap=-2047 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      common/z1in/leg(1:maxleg),nleg,nloop,incom
      common/z3in/momep(0:maxleg),momel(0:maxleg)
      common/z5in/dunit,munit,ounit,sunit,funit
      common/z2g/dis,dsym
      common/z3g/nv(0:0),dpntro(0:0),vparto(0:0),vdeg(0:0),nvert
      common/z4g/n,nli
      common/z6g/p1(1:maxleg),invp1(1:maxleg)
      common/z7g/lmap(1:maxn,1:maxdeg),vmap(1:maxn,1:maxdeg),
     :pmap(1:maxn,1:maxdeg),vlis(1:maxn),invlis(1:maxn)
      common/z8g/degree(1:maxn),xn(1:maxn)
      common/z11g/nphi,nblok,nprop,npprop
      common/z16g/rdeg(1:maxn),amap(1:maxn,1:maxdeg)
      common/z18g/eg(1:maxn,1:maxn),flow(1:maxli,0:maxleg+maxrho)
      common/z19g/vfo(1:maxn)
      character*(srec) auxlin
      common/z22g/auxlin
      common/z24g/iogp(1:4)
      common/z25g/ex(1:maxli),ey(1:maxli),ovm(1:maxn,1:maxdeg)
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z33g/namep(0:0),namel(0:0)
      common/z34g/link(0:0),antiq(0:0),blok(0:0)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      common/z38g/gmkp(0:0),gmkl(0:0),gmkd(0:0),gmko(0:0),
     :gmkvp(0:0),gmkvl(0:0)
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      common/z47g/ndiagp,ndiagl,hhp,hhl,doffp,doffl,noffp,noffl,wsint
      common/z51g/aplus,aminus,azero,anine,aeq
      common/z52g/acf0(0:127),acf1(0:127),lfeed,aspace,dquote,squote
      ig=iogp(2)
      klin=0
      lupt=0
   10 continue
      if((ig.lt.iogp(2)).or.(ig.ge.iogp(3)))then
      goto 80
      endif
      k=stib(ig)
      if(k.gt.0)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(k)
      ig=ig+1
      goto 10
      endif
      ii=stib(ig+1)
      if(k.gt.-3)then
      if(k.eq.-1)then
      lupt=lupty(ii)
      if(lupt.gt.0)then
      if(stib(ig+4).eq.0)then
      if(lupt.eq.3)then
      if(ii.eq.13)then
      stib(ig+4)=incom
      stib(ig+5)=0
      elseif(ii.eq.14)then
      stib(ig+4)=nleg
      stib(ig+5)=incom
      endif
      else
      if(lupt.eq.4)then
      stib(ig+4)=nli-nleg
      elseif(lupt.eq.5)then
      stib(ig+4)=n-nleg
      elseif(lupt.eq.6)then
      stib(ig+4)=degree(nleg+lupi)
      endif
      stib(ig+5)=0
      endif
      endif
      stib(ig+5)=stib(ig+5)+1
      if(stib(ig+4).lt.stib(ig+5))then
      stib(ig+4)=0
      stib(ig+5)=0
      if(lupt.lt.6)then
      lupt=0
      else
      lupt=5
      endif
      ig=ig+1
      else
      if(lupt.lt.6)then
      lupi=stib(ig+5)
      else
      lupj=stib(ig+5)
      endif
      endif
      elseif(lupt.eq.0)then
      goto 80
      endif
      elseif(k.eq.-2)then
      if(ii.eq.21)then
      if(klin.gt.0)then
      klin=klin-1
      else
      goto 80
      endif
      else
      goto 80
      endif
      else
      goto 80
      endif
      elseif(k.gt.-5)then
      if(k.eq.-3)then
      if(ii.gt.70.and.ii.lt.80)then
      if(ii.lt.74)then
      if(ii.eq.72)then
      jj=dsym
      if(dsym.gt.1)then
      klin=klin+2
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip-1:ip)='1/'
      endif
      elseif(ii.eq.73)then
      jj=dsym
      endif
      elseif(ii.lt.77)then
      if(ii.eq.74)then
      jj=nli-nleg
      elseif(ii.eq.75)then
      jj=nleg
      elseif(ii.eq.76)then
      jj=nloop
      endif
      else
      if(ii.eq.77)then
      jj=n-nleg
      elseif(ii.eq.78)then
      jj=incom
      elseif(ii.eq.79)then
      jj=nleg-incom
      endif
      endif
      ip=stcbs(1)+klin
      if(ii.ne.71)then
      call dkar(jj,ip,jk)
      klin=klin+jk
      else
      klin=klin+noffl
      call vaocb(klin)
      stcb(ip+1:ip+noffl)=stcb(noffp+1:noffp+noffl)
      endif
      else
      if(ii.eq.61)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(44-dis)
      elseif(ii.eq.62)then
      if(dis.lt.0)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(aminus)
      endif
      else
      goto 80
      endif
      endif
      elseif(k.eq.-4)then
      if(lupt.eq.4)then
      if(ii.lt.40)then
      if(ii.eq.31.or.ii.eq.33)then
      i=lupi
      j=pmap(ex(i),ey(i))
      if(ii.eq.33)then
      j=stib(link(0)+j)
      endif
      il=stib(namel(0)+j)
      klin=klin+il
      call vaocb(klin)
      ia=stib(namep(0)+j)
      ip=stcbs(1)+klin
      stcb(ip-il+1:ip)=stcb(ia:ia-1+il)
      elseif(ii.eq.32.or.ii.eq.34)then
      i=lupi
      if(ii.eq.32)then
      ik=1
      else
      ik=0
      endif
      j=ey(i)
      i=ex(i)
      k0=klin
      other=0
      if(vmap(i,j).eq.i)then
      if(mod(j-rdeg(i)+ik,2).ne.0)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(aminus)
      endif
      klin=klin+momel(0)
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip-momel(0)+1:ip)=
     :stcb(momep(0):momep(0)-1+momel(0))
      ij=eg(i,i)+(j-1-rdeg(i))/2
      do 372 k=nleg+1,nleg+nloop
      if(flow(ij,k).ne.0)then
      ip=stcbs(1)+klin
      call dkar(k-nleg,ip,jk)
      klin=klin+jk
      goto 70
      endif
  372 continue
      else
      if(vmap(i,j).lt.i)then
      ij=-1
      else
      ij=1
      endif
      if(ik.eq.0)then
      ij=-ij
      endif
      do 371 k=nleg+1,nleg+nloop
      jj=flow(amap(i,j),k)
      if(jj.ne.0)then
      if(other.eq.1.or.ij.eq.jj)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(44+ij*jj)
      endif
      klin=klin+momel(0)
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip-momel(0)+1:ip)=
     :stcb(momep(0):momep(0)-1+momel(0))
      call dkar(k-nleg,ip,jk)
      klin=klin+jk
      other=1
      endif
  371 continue
      do 363 k=1,nleg
      jj=flow(amap(i,j),invp1(k))
      if(jj.ne.0)then
      if(k.gt.incom)then
      jj=-jj
      endif
      if(other.eq.1.or.ij.eq.jj)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(44+ij*jj)
      endif
      klin=klin+momel(k)
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip-momel(k)+1:ip)=
     :stcb(momep(k):momep(k)-1+momel(k))
      other=1
      endif
  363 continue
      if(k0.eq.klin)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(azero)
      endif
      endif
      elseif(ii.eq.35)then
      i=lupi
      j=pmap(ex(i),ey(i))
      klin=klin+1
      call vaocb(klin)
      ia=stib(antiq(0)+j)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(43+ia+ia)
      else
      goto 80
      endif
      elseif(ii.lt.50)then
      if(ii.lt.44)then
      if(ii.eq.40)then
      jj=3
      elseif(ii.eq.41)then
      jj=lupi
      elseif(ii.eq.42)then
      jj=2*lupi-1
      elseif(ii.eq.43)then
      i=lupi
      j=ey(i)
      i=ex(i)
      jj=0
      do 40 k=1,degree(i)
      if(ovm(i,k).eq.j)then
      jj=k
      endif
   40 continue
      else
      goto 80
      endif
      elseif(ii.lt.47)then
      if(ii.eq.44)then
      i=lupi
      jj=ex(i)-nleg
      elseif(ii.eq.45)then
      jj=2*lupi
      elseif(ii.eq.46)then
      i=lupi
      j=lmap(ex(i),ey(i))
      i=vmap(ex(i),ey(i))
      jj=0
      do 50 k=1,degree(i)
      if(ovm(i,k).eq.j)then
      jj=k
      endif
   50 continue
      else
      goto 80
      endif
      elseif(ii.lt.50)then
      if(ii.eq.47)then
      i=lupi
      jj=vmap(ex(i),ey(i))-nleg
      elseif(ii.eq.48)then
      i=lupi
      jj=degree(ex(i))
      elseif(ii.eq.49)then
      i=lupi
      jj=degree(vmap(ex(i),ey(i)))
      else
      goto 80
      endif
      else
      goto 80
      endif
      ip=stcbs(1)+klin
      call dkar(jj,ip,jk)
      klin=klin+jk
      else
      goto 80
      endif
      elseif(lupt.eq.3)then
      if(ii.lt.40)then
      if(ii.eq.31.or.ii.eq.33)then
      i=lupi
      j=stib(link(0)+leg(i))
      if(ii.eq.33)then
      j=stib(link(0)+j)
      endif
      if(lupi.gt.incom)then
      j=stib(link(0)+j)
      endif
      il=stib(namel(0)+j)
      klin=klin+il
      call vaocb(klin)
      ia=stib(namep(0)+j)
      ip=stcbs(1)+klin
      stcb(ip-il+1:ip)=stcb(ia:ia-1+il)
      elseif(ii.eq.32.or.ii.eq.34)then
      i=lupi
      if(ii.eq.34)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(aminus)
      endif
      klin=klin+momel(i)
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip-momel(i)+1:ip)=
     :stcb(momep(i):momep(i)-1+momel(i))
      elseif(ii.eq.35)then
      i=lupi
      j=leg(i)
      klin=klin+1
      call vaocb(klin)
      ia=stib(antiq(0)+j)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(43+ia+ia)
      else
      goto 80
      endif
      elseif(ii.lt.50)then
      if(ii.eq.40)then
      i=lupi
      if(i.le.incom)then
      jj=1
      else
      jj=2
      endif
      elseif(ii.eq.42)then
      i=lupi
      if(i.le.incom)then
      jj=1-2*i
      else
      jj=2*(incom-i)
      endif
      elseif(ii.eq.43)then
      i=lupi
      j=lmap(invp1(i),1)
      i=vmap(invp1(i),1)
      jj=0
      do 44 k=1,degree(i)
      if(ovm(i,k).eq.j)then
      jj=k
      endif
   44 continue
      elseif(ii.eq.44)then
      i=lupi
      jj=vmap(invp1(i),1)-nleg
      elseif(ii.eq.48)then
      i=lupi
      j=lmap(invp1(i),1)
      jj=degree(vmap(invp1(i),1))
      else
      goto 80
      endif
      ip=stcbs(1)+klin
      call dkar(jj,ip,jk)
      klin=klin+jk
      elseif(ii.lt.60)then
      if(ii.eq.51)then
      i=lupi
      elseif(ii.eq.52)then
      i=lupi
      elseif(ii.eq.53)then
      i=lupi-incom
      else
      goto 80
      endif
      ip=stcbs(1)+klin
      call dkar(i,ip,jk)
      klin=klin+jk
      else
      goto 80
      endif
      elseif(lupt.eq.5.or.lupt.eq.6)then
      if(ii.lt.40)then
      if(ii.eq.31.or.ii.eq.33)then
      i=nleg+lupi
      j=lupj
      j=pmap(i,ovm(i,j))
      if(ii.eq.33)then
      j=stib(link(0)+j)
      endif
      il=stib(namel(0)+j)
      klin=klin+il
      call vaocb(klin)
      ia=stib(namep(0)+j)
      ip=stcbs(1)+klin
      stcb(ip-il+1:ip)=stcb(ia:ia-1+il)
      elseif(ii.eq.32.or.ii.eq.34)then
      i=nleg+lupi
      j=lupj
      j=ovm(i,j)
      ij=vmap(i,j)
      if(ii.eq.32)then
      kk=0
      else
      kk=1
      endif
      k0=klin
      other=0
      if(ij.le.nleg)then
      ij=p1(ij)
      if(ij.gt.incom)then
      kk=1-kk
      endif
      if(kk.eq.1)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(aminus)
      endif
      klin=klin+momel(ij)
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip-momel(ij)+1:ip)=
     :stcb(momep(ij):momep(ij)-1+momel(ij))
      elseif(vmap(i,j).eq.i)then
      if(mod(j-rdeg(i)+kk,2).eq.0)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(aminus)
      endif
      klin=klin+momel(0)
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip-momel(0)+1:ip)=
     :stcb(momep(0):momep(0)-1+momel(0))
      ij=eg(i,i)+(j-1-rdeg(i))/2
      do 172 k=nleg+1,nleg+nloop
      if(flow(ij,k).ne.0)then
      ip=stcbs(1)+klin
      call dkar(k-nleg,ip,jk)
      klin=klin+jk
      goto 70
      endif
  172 continue
      else
      ij=1
      if(vmap(i,j).lt.i)then
      ij=-ij
      endif
      if(kk.eq.1)then
      ij=-ij
      endif
      do 171 k=nleg+1,nleg+nloop
      jj=flow(amap(i,j),k)
      if(jj.ne.0)then
      if(other.eq.1.or.ij.eq.jj)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(44+ij*jj)
      endif
      klin=klin+momel(0)
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip-momel(0)+1:ip)=
     :stcb(momep(0):momep(0)-1+momel(0))
      call dkar(k-nleg,ip,jk)
      klin=klin+jk
      other=1
      endif
  171 continue
      do 163 k=1,nleg
      jj=flow(amap(i,j),invp1(k))
      if(jj.ne.0)then
      if(k.gt.incom)then
      jj=-jj
      endif
      if(other.eq.1.or.ij.eq.jj)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(44+ij*jj)
      endif
      klin=klin+momel(k)
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip-momel(k)+1:ip)=
     :stcb(momep(k):momep(k)-1+momel(k))
      other=1
      endif
  163 continue
      if(k0.eq.klin)then
      klin=klin+1
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(azero)
      endif
      endif
      elseif(ii.eq.35)then
      i=nleg+lupi
      j=lupj
      j=pmap(i,ovm(i,j))
      klin=klin+1
      call vaocb(klin)
      ia=stib(antiq(0)+j)
      ip=stcbs(1)+klin
      stcb(ip:ip)=char(43+ia+ia)
      endif
      elseif(ii.lt.50)then
      if(ii.eq.41)then
      i=nleg+lupi
      j=lupj
      k=stib(stib(vparto(0)+stib(vfo(i)))+j)
      j=ovm(i,j)
      ij=vmap(i,j)
      if(ij.le.nleg)then
      ij=p1(ij)
      if(ij.le.incom)then
      jj=1-2*ij
      else
      jj=2*(incom-ij)
      endif
      else
      jj=amap(i,j)-nleg
      endif
      elseif(ii.lt.46)then
      if(ii.eq.40)then
      i=nleg+lupi
      j=lupj
      k=stib(stib(vparto(0)+stib(vfo(i)))+j)
      j=ovm(i,j)
      ij=vmap(i,j)
      if(ij.le.nleg)then
      if(p1(ij).le.incom)then
      jj=1
      else
      jj=2
      endif
      else
      jj=3
      endif
      elseif(ii.eq.42.or.ii.eq.45)then
      i=nleg+lupi
      j=lupj
      k=stib(stib(vparto(0)+stib(vfo(i)))+j)
      j=ovm(i,j)
      ij=vmap(i,j)
      if(ij.le.nleg)then
      ij=p1(ij)
      if(ij.le.incom)then
      jj=1-2*ij
      else
      jj=2*(incom-ij)
      endif
      elseif(k.lt.stib(link(0)+k))then
      jj=2*(amap(i,j)-nleg)-1
      elseif(k.gt.stib(link(0)+k))then
      jj=2*(amap(i,j)-nleg)
      elseif(i.lt.ij)then
      jj=2*(amap(i,j)-nleg)-1
      elseif(i.gt.ij)then
      jj=2*(amap(i,j)-nleg)
      elseif(mod(j-rdeg(i),2).ne.0)then
      jj=2*(amap(i,j)-nleg)-1
      else
      jj=2*(amap(i,j)-nleg)
      endif
      if(ii.eq.45)then
      if(jj.gt.0)then
      jj=jj-1+2*mod(jj,2)
      else
      jj=0
      endif
      endif
      elseif(ii.eq.43)then
      jj=lupj
      elseif(ii.eq.44)then
      jj=lupi
      else
      goto 80
      endif
      elseif(ii.lt.50)then
      if(ii.eq.46)then
      i=nleg+lupi
      j=lupj
      ij=vmap(i,ovm(i,j))
      j=lmap(i,ovm(i,j))
      jj=0
      if(ij.gt.nleg)then
      do 74 k=1,degree(ij)
      if(ovm(ij,k).eq.j)then
      jj=k
      endif
   74 continue
      endif
      elseif(ii.eq.47)then
      i=nleg+lupi
      j=lupj
      jj=vmap(i,ovm(i,j))-nleg
      if(jj.lt.0)then
      jj=0
      endif
      elseif(ii.eq.48)then
      jj=degree(nleg+lupi)
      elseif(ii.eq.49)then
      i=nleg+lupi
      j=lupj
      jj=vmap(i,ovm(i,j))
      if(jj.le.nleg)then
      jj=0
      else
      jj=degree(jj)
      endif
      else
      goto 80
      endif
      else
      goto 80
      endif
      ip=stcbs(1)+klin
      call dkar(jj,ip,jk)
      klin=klin+jk
      elseif(ii.lt.60)then
      else
      goto 80
      endif
      else
      goto 80
      endif
      else
      goto 80
      endif
      elseif(k.eq.-6)then
      if(ii.le.0.or.ii.gt.nudk)then
      goto 80
      endif
      if(stib(udkt(0)+ii).eq.1)then
      ij=stib(udki(0)+ii)
      if(ij.le.0.or.ij.gt.ngmk)then
      goto 80
      elseif(stib(gmkd(0)+ij).ne.1)then
      goto 80
      endif
      ij=stib(gmko(0)+ij)+1
      i=stib(gmkvp(0)+ij)
      j=stib(gmkvl(0)+ij)
      if(j.gt.0)then
      klin=klin+j
      call vaocb(klin)
      ip=stcbs(1)+klin
      stcb(ip-j+1:ip)=stcb(i:i-1+j)
      endif
      if(j.lt.0.or.stib(ig+3).ne.0)then
      goto 80
      endif
      elseif(stib(udkt(0)+ii).eq.2)then
      if(lupt.eq.3)then
      i=lupi
      j=leg(i)
      if(stib(ig+3).eq.0)then
      j=stib(link(0)+j)
      endif
      if(lupi.gt.incom)then
      j=stib(link(0)+j)
      endif
      elseif(lupt.eq.4)then
      i=lupi
      j=pmap(ex(i),ey(i))
      if(stib(ig+3).ne.0)then
      j=stib(link(0)+j)
      endif
      elseif(lupt.eq.5)then
      goto 80
      elseif(lupt.eq.6)then
      i=nleg+lupi
      j=lupj
      j=pmap(i,ovm(i,j))
      if(stib(ig+3).ne.0)then
      j=stib(link(0)+j)
      endif
      else
      goto 80
      endif
      ij=stib(udki(0)+ii)
      if(ij.le.0.or.ij.gt.npmk)then
      goto 80
      elseif(j.le.0.or.j.gt.nphi)then
      goto 80
      endif
      il=stib(stib(pmkvlp(0)+ij)+j)
      if(il.gt.0)then
      klin=klin+il
      call vaocb(klin)
      ia=stib(stib(pmkvpp(0)+ij)+j)
      ip=stcbs(1)+klin
      stcb(ip-il+1:ip)=stcb(ia:ia-1+il)
      endif
      if(il.lt.0)then
      goto 80
      endif
      elseif(stib(udkt(0)+ii).eq.3)then
      if(lupt.eq.3)then
      i=lupi
      jj=vmap(invp1(i),1)
      elseif(lupt.eq.4)then
      i=lupi
      if(stib(ig+3).eq.0)then
      jj=ex(i)
      else
      jj=vmap(ex(i),ey(i))
      endif
      elseif(lupt.eq.5)then
      jj=nleg+lupi
      elseif(lupt.eq.6)then
      jj=nleg+lupi
      else
      goto 80
      endif
      j=stib(vfo(jj))
      ij=stib(udki(0)+ii)
      if(ij.le.0.or.ij.gt.nvmk)then
      goto 80
      elseif(j.le.0.or.j.gt.nvert)then
      goto 80
      endif
      il=stib(stib(vmkvlp(0)+ij)+j)
      if(il.gt.0)then
      klin=klin+il
      call vaocb(klin)
      ia=stib(stib(vmkvpp(0)+ij)+j)
      ip=stcbs(1)+klin
      stcb(ip-il+1:ip)=stcb(ia:ia-1+il)
      endif
      if(il.lt.0)then
      goto 80
      endif
      else
      goto 80
      endif
      elseif(k.eq.eoa)then
      goto 20
      else
      goto 80
      endif
   70 continue
      ig=ig+stib(ig+2)
      goto 10
   20 continue
      j1=stcbs(1)+1
      do i1=stcbs(1)+1,stcbs(1)+klin
      if(ichar(stcb(i1:i1)).eq.lfeed)then
      ios=1
      if(i1.gt.j1)then
      write(unit=ounit,fmt=404,iostat=ios)stcb(j1:i1-1)
      else
      write(unit=ounit,fmt=404,iostat=ios)
      endif
      if(ios.ne.0)then
      auxlin(1:srec)='run-time error while writing to'
      call messag(1,0,0,5)
      endif
      j1=i1+1
      endif
      enddo
  404 format(a)
      goto 90
   80 continue
      auxlin(1:srec)='run-time error while processing <diagram>'
      call messag(1,0,0,0)
   90 continue
      return
      end
      subroutine style
      implicit integer(a-z)
      save
      parameter ( maxtak=2 )
      parameter ( maxleg=11, maxrho=7, maxi=11, maxdeg=6 )
      parameter ( srec=81, ssrec=62 )
      parameter ( eoa=-127, nap=-2047 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( sxbuff=2040 )
      common/z5in/dunit,munit,ounit,sunit,funit
      common/z7in/lsfile,sfilea,sfileb
      character*(srec) auxlin
      common/z22g/auxlin
      common/z23g/tak(1:maxtak),ks
      common/z24g/iogp(1:4)
      common/z26g/kes(0:0),kle(0:0),pstke(0:0),wstke(0:0)
      character*(sxbuff) stxb
      common/z27g/stxb
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      common/z52g/acf0(0:127),acf1(0:127),lfeed,aspace,dquote,squote
      integer lg(srec)
      ks=0
      level=0
      bevel=0
      do i1=1,srec
      lg(i1)=0
      enddo
      nudk=0
      call aoib(2)
      udkp(1)=stibs(1)-1
      stib(udkp(1))=eoa
      stib(stibs(1))=eoa
      sunit=4
      call iopen(sfilea,lsfile,sunit,0)
      nlin=0
   70 continue
      call qread(sunit,4,nlin,slin)
      if(slin.eq.-1)then
      goto 444
      endif
      if(slin.gt.0)then
      stxb(1:slin)=stcb(stcbs(1)+1:stcbs(1)+slin)
      endif
      if(abs(level-2).eq.2)then
      ii=0
      if(slin.eq.0)then
      ii=1
      elseif(stxb(1:1).eq.'#')then
      ii=1
      elseif(stxb(1:1).eq.'%')then
      ii=1
      elseif(stxb(1:1).eq.'*')then
      ii=1
      endif
      if(ii.ne.0)then
      goto 70
      endif
      endif
      if(level.eq.4)then
      if(slin.gt.0)then
      goto 300
      endif
      endif
      iw=slin
      if(level.eq.0)then
      if(stxb(1:1).ne.'<')then
      goto 300
      endif
      endif
      j1=0
      j2=0
      k1=0
      k2=0
      nlg=0
      nlg2=0
      do 27 i=1,iw
      if(j2.eq.0.and.stxb(i:i).eq.'<')then
      j1=j1+1
      if(j1.eq.2)then
      call aoib(1)
      stib(stibs(1))=ichar(stxb(i:i))
      j1=0
      elseif(stxb(i+1:i+1).ne.'<')then
      if(j2.ne.0)then
      goto 300
      endif
      nlg=nlg+1
      nlg2=1-nlg2
      if(nlg2.ne.1)then
      goto 300
      endif
      lg(nlg)=i
      endif
      elseif(j1.eq.0.and.stxb(i:i).eq.'[')then
      j2=j2+1
      if(j2.eq.2)then
      call aoib(1)
      stib(stibs(1))=ichar(stxb(i:i))
      j2=0
      elseif(stxb(i+1:i+1).ne.'[')then
      if(j1.ne.0)then
      goto 300
      endif
      nlg=nlg+1
      nlg2=1-nlg2
      if(nlg2.ne.1)then
      goto 300
      endif
      lg(nlg)=i
      endif
      elseif(j2.eq.0.and.stxb(i:i).eq.'>')then
      if(j1.eq.0)then
      k1=k1+1
      if(k1.eq.2)then
      call aoib(1)
      stib(stibs(1))=ichar(stxb(i:i))
      k1=0
      elseif(stxb(i+1:i+1).ne.'>')then
      goto 300
      endif
      else
      j1=0
      if(j2.ne.0)then
      goto 300
      endif
      nlg=nlg+1
      nlg2=1-nlg2
      if(nlg2.ne.0)then
      goto 300
      endif
      lg(nlg)=i
      ii=lg(nlg)-lg(nlg-1)-1
      if(ii.le.0)then
      goto 300
      endif
      ii=styki(stxb,lg(nlg-1)+1,lg(nlg)-1,level)
      if(ii.eq.0)then
      goto 300
      endif
      if(ii.gt.0.and.ii.lt.5)then
      if(nlg.ne.2.or.lg(1).ne.1.or.lg(2).ne.iw)then
      goto 300
      endif
      if(ks.ne.0)then
      goto 300
      endif
      level=ii
      else
      call ktoc(ii,nlin)
      endif
      endif
      elseif(j1.eq.0.and.stxb(i:i).eq.']')then
      if(j2.eq.0)then
      k2=k2+1
      if(k2.eq.2)then
      call aoib(1)
      stib(stibs(1))=ichar(stxb(i:i))
      k2=0
      elseif(stxb(i+1:i+1).ne.']')then
      goto 300
      endif
      else
      j2=0
      if(j1.ne.0)then
      goto 300
      endif
      nlg=nlg+1
      nlg2=1-nlg2
      if(nlg2.ne.0)then
      goto 300
      endif
      lg(nlg)=i
      lw=lg(nlg)-lg(nlg-1)-1
      if(lw.le.0)then
      goto 300
      endif
      ii=udstyk(stxb,lg(nlg-1)+1,lg(nlg)-1)
      if(ii.eq.0)then
      goto 300
      endif
      endif
      elseif(j1.eq.0.and.j2.eq.0)then
      call aoib(1)
      stib(stibs(1))=ichar(stxb(i:i))
      endif
   27 continue
      if(nlg2.ne.0)then
      goto 300
      endif
      if(level.gt.1.or.bevel.eq.level)then
      call aoib(1)
      endif
      if(bevel.eq.level)then
      stib(stibs(1))=lfeed
      else
      bevel=level
      if(level.gt.1)then
      stib(stibs(1))=eoa
      endif
      iogp(level)=stibs(1)+1
      endif
      goto 70
  444 continue
      if(level.lt.4)then
      auxlin(1:srec)='is incomplete'
      call messag(1,0,0,-4)
      endif
      if(ks.eq.0)then
      call rsfki
      goto 90
      endif
  300 continue
      auxlin(1:srec)='wrong syntax,'
      call messag(1,nlin,nlin,4)
   90 continue
      return
      end
      subroutine ktoc(ii,nlin)
      implicit integer(a-z)
      save
      parameter ( maxtak=2 )
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      character*(srec) auxlin
      common/z22g/auxlin
      common/z23g/tak(1:maxtak),ks
      common/z30g/stib(1:sibuff)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z52g/acf0(0:127),acf1(0:127),lfeed,aspace,dquote,squote
      integer astak(1:maxtak)
      if((ii.lt.10).or.(ii.gt.82))then
      auxlin(1:srec)='ktoc_1'
      call messag(1,0,0,0)
      endif
      if(ii.le.30)then
      if(ii.lt.18)then
      if(ks.ge.maxtak)then
      auxlin(1:srec)='too many nested loops,'
      call messag(1,nlin,nlin,4)
      endif
      call aoib(6)
      stib(stibs(1)-5)=-1
      stib(stibs(1)-4)=ii
      stib(stibs(1)-3)=6
      stib(stibs(1)-2)=0
      stib(stibs(1)-1)=0
      stib(stibs(1))=0
      ks=ks+1
      tak(ks)=ii
      astak(ks)=stibs(1)-5
      if(ii.eq.11.or.ii.eq.12)then
      if(ii.eq.11)then
      if(ks.ne.1)then
      goto 111
      endif
      else
      if(ks.ne.2)then
      goto 111
      elseif(tak(ks-1).ne.11)then
      goto 111
      endif
      endif
      elseif(ii.lt.17)then
      if(ks.ne.1)then
      goto 111
      endif
      elseif(ii.eq.17)then
      if(ks.ne.2)then
      goto 111
      elseif(tak(ks-1).ne.16)then
      goto 111
      endif
      endif
      elseif(ii.eq.19)then
      if(ks.le.0)then
      goto 111
      endif
      call aoib(3)
      stib(stibs(1)-2)=-1
      stib(stibs(1)-1)=ii
      jj=stibs(1)-astak(ks)
      stib(stibs(1))=2-jj
      stib(astak(ks)+3)=jj
      ks=ks-1
      elseif(ii.eq.21)then
      call aoib(3)
      stib(stibs(1)-2)=-2
      stib(stibs(1)-1)=ii
      stib(stibs(1))=3
      elseif(ii.eq.22)then
      call aoib(1)
      stib(stibs(1))=lfeed
      else
      goto 111
      endif
      endif
      if(ii.gt.30)then
      if(ii.lt.60.or.ii.eq.82)then
      call aoib(3)
      stib(stibs(1)-2)=-4
      stib(stibs(1)-1)=ii
      stib(stibs(1))=3
      elseif(ii.lt.90)then
      call aoib(3)
      stib(stibs(1)-2)=-3
      stib(stibs(1)-1)=ii
      stib(stibs(1))=3
      endif
      endif
      if(ii.eq.82)then
      if(ks.ne.1.and.ks.ne.2)then
      goto 111
      endif
      endif
      if(ii.lt.31.or.ii.gt.60)then
      goto 90
      endif
      j=0
      if(ks.gt.0)then
      j=tak(ks)
      endif
      if(styk(ii,j).eq.0)then
      goto 111
      endif
      goto 90
  111 continue
      auxlin(1:srec)='wrong syntax,'
      call messag(1,nlin,nlin,4)
   90 continue
      return
      end
      integer function styk(k,pl)
      implicit integer(a-z)
      save
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      common/z13g/kloo(15:32,1:4),loopt(11:19)
      common/z26g/kes(0:0),kle(0:0),pstke(0:0),wstke(0:0)
      common/z30g/stib(1:sibuff)
      if(k.lt.stib(kes(0)+15))then
      goto 200
      elseif(k.gt.stib(kes(0)+32))then
      goto 200
      elseif(lupty(pl).eq.0)then
      goto 200
      endif
      do 110 i=15,32
      if(stib(kes(0)+i).eq.k)then
      styk=kloo(i,lupty(pl)-2)
      goto 300
      endif
  110 continue
  200 continue
      styk=0
  300 continue
      if(k.eq.52)then
      if(pl.ne.13)then
      styk=0
      endif
      elseif(k.eq.53.or.k.eq.51)then
      if(pl.ne.14)then
      styk=0
      endif
      endif
      return
      end
      integer function styki(s,s1,s2,level)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      character*(srec) s
      common/z12g/jflag(1:11),mflag(1:24)
      character*(srec) auxlin
      common/z22g/auxlin
      common/z26g/kes(0:0),kle(0:0),pstke(0:0),wstke(0:0)
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      if(s1.gt.s2)then
      auxlin(1:srec)='styki_1'
      call messag(1,0,0,0)
      endif
      styki=0
      call mstr0(s,s1,s2,pstke(0),wstke(0),jk)
      if(jk.ne.0)then
      styki=stib(kes(0)+jk)
      if(styki.le.4)then
      if(styki.ne.level+1)then
      styki=0
      endif
      elseif(level.le.3)then
      k=1
      do 10 i=2,level
      k=k*2
   10 continue
      if(mod(stib(kle(0)+jk),2*k).lt.k)then
      styki=0
      endif
      endif
      endif
      if(styki.eq.32.or.styki.eq.34)then
      jflag(10)=1
      endif
      return
      end
      integer function udstyk(s,s1,s2)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( eoa=-127, nap=-2047 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( maxtak=2 )
      common/z23g/tak(1:maxtak),ks
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      common/z45g/qdatp,qdatl,qvp,qvl,dprefp,dprefl
      common/z52g/acf0(0:127),acf1(0:127),lfeed,aspace,dquote,squote
      character*(srec) s
      kl=0
      if(ks.gt.0)then
      kl=lupty(tak(ks))
      endif
      udstyk=0
      if(kl.lt.3.or.kl.gt.6)then
      if(kl.ne.0)then
      goto 90
      endif
      kl=2
      endif
      kd=0
      j1=s1-1+dprefl
      if(s2.ge.j1)then
      if(s(s1:j1).eq.stcb(dprefp:dprefp-1+dprefl))then
      if(s2.gt.j1)then
      kd=1
      else
      goto 90
      endif
      endif
      endif
      j1=s1+kd*dprefl
      jj=s2-j1+1
      if(stdw(s,j1,s2).eq.0)then
      goto 90
      endif
      if(nudk.gt.0)then
      call mstr1(s,j1,s2,udkp(1),1,2,udstyk)
      endif
      j2=stibs(1)
      kstep=4
      call aoib(kstep)
      stib(j2+1)=-6
      if(udstyk.eq.0)then
      nudk=nudk+1
      stib(j2+2)=nudk
      else
      stib(j2+2)=udstyk
      endif
      stib(j2+3)=kstep
      stib(j2+4)=kd
      ii=udkp(1)
      if(udstyk.ne.0)then
      do i1=1,udstyk
      ii=stib(ii)
      enddo
      j2=ii+1+kl
      if(stib(j2).eq.0)then
      stib(j2)=1+kd
      elseif(stib(j2).eq.2-kd)then
      stib(j2)=3
      endif
      goto 90
      endif
      udstyk=stib(j2+2)
      do i1=2,nudk
      ii=stib(ii)
      enddo
      j2=stibs(1)
      kstep=8
      call aoib(kstep)
      stib(j2-1)=stib(j2-1)+kstep
      stib(ii)=j2+1
      stib(j2+1)=eoa
      call aocb(jj+1)
      stib(j2+2)=stcbs(1)-jj
      stcb(stib(j2+2):stcbs(1))=s(j1:s2)//char(lfeed)
      stib(j2+3)=jj
      stib(j2+4)=0
      stib(j2+5)=0
      stib(j2+6)=0
      stib(j2+7)=0
      stib(j2+8)=0
      stib(j2+2+kl)=kd+1
   90 continue
      return
      end
      subroutine stpa(s,ia,ib,ds)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      common/z21g/punct1(0:0)
      character*(srec) auxlin
      common/z22g/auxlin
      common/z30g/stib(1:sibuff)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z52g/acf0(0:127),acf1(0:127),lfeed,aspace,dquote,squote
      character*(*) s
      if(ia.le.0.or.ib.lt.ia.or.ds.lt.0)then
      auxlin(1:srec)='stpa_1'
      call messag(1,0,0,0)
      endif
      j1=0
      j2=0
      it=0
      plic=0
      sts=ds+1
      call vaoib(sts)
      st0=stibs(1)+sts
      s0=st0
      stib(s0)=0
      do i1=ia,ib+1
      if(i1.le.ib)then
      ic=ichar(s(i1:i1))
      else
      ic=lfeed
      endif
      if(ic.eq.squote)then
      plic=1-plic
      if(plic.eq.1)then
      if(it.eq.0)then
      it=2
      j1=i1
      endif
      endif
      j2=i1
      elseif(plic.eq.1)then
      if(ic.eq.lfeed)then
      stib(s0)=-1
      endif
      j2=i1
      elseif(stib(punct1(0)+ic).ne.0)then
      if(it.eq.2)then
      if(stds(s,j1,j2,0,1).ge.0)then
      it=2
      elseif(stdw(s,j1,j2).ne.0)then
      it=3
      elseif(stdz(s,j1,j2).ne.0)then
      it=4
      elseif(stdq(s,j1,j2).ne.0)then
      it=5
      else
      stib(s0)=-1
      goto 70
      endif
      sts=sts+3
      call vaoib(sts)
      st0=st0+3
      stib(st0-2)=it
      stib(st0-1)=j1
      stib(st0)=j2
      stib(s0)=stib(s0)+1
      j1=0
      it=0
      endif
      sts=sts+3
      call vaoib(sts)
      st0=st0+3
      stib(st0-2)=1
      stib(st0-1)=stib(punct1(0)+ic)
      stib(st0)=i1
      if(stib(s0).ge.0)then
      stib(s0)=stib(s0)+1
      endif
      elseif((ic.eq.lfeed).or.(ic.eq.aspace))then
      if(it.eq.2)then
      if(stds(s,j1,j2,0,1).ge.0)then
      it=2
      elseif(stdw(s,j1,j2).ne.0)then
      it=3
      elseif(stdz(s,j1,j2).ne.0)then
      it=4
      elseif(stdq(s,j1,j2).ne.0)then
      it=5
      else
      stib(s0)=-1
      goto 70
      endif
      sts=sts+3
      call vaoib(sts)
      st0=st0+3
      stib(st0-2)=it
      stib(st0-1)=j1
      stib(st0)=j2
      stib(s0)=stib(s0)+1
      j1=0
      it=0
      endif
      else
      if(it.eq.0)then
      it=2
      j1=i1
      elseif(it.ne.2)then
      stib(s0)=-1
      endif
      j2=i1
      endif
   70 continue
      if(stib(s0).eq.-1)then
      goto 90
      endif
      enddo
      if(plic.ne.0)then
      stib(s0)=-1
      endif
   90 continue
      return
      end
      subroutine hsort(x,ia,ib)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      character*(srec) auxlin
      common/z22g/auxlin
      integer x(*)
      if((ia.le.0).or.(ia.gt.ib))then
      auxlin(1:srec)='hsort_1'
      call messag(1,0,0,0)
      elseif(ia.eq.ib)then
      goto 90
      endif
      i0=ia-1
      hr=ib-i0
      hl=1+hr/2
   20 continue
      if(hl.gt.1)then
      hl=hl-1
      xtmp=x(i0+hl)
      else
      xtmp=x(i0+hr)
      x(i0+hr)=x(i0+1)
      hr=hr-1
      if(hr.eq.1)then
      x(i0+1)=xtmp
      goto 90
      endif
      endif
      hj=hl
   40 continue
      hi=hj
      hj=hi+hi
      if(hj.le.hr)then
      if(hj.lt.hr)then
      if(x(i0+hj).lt.x(i0+hj+1))then
      hj=hj+1
      endif
      endif
      if(xtmp.lt.x(i0+hj))then
      x(i0+hi)=x(i0+hj)
      goto 40
      endif
      endif
      x(i0+hi)=xtmp
      goto 20
   90 continue
      return
      end
      subroutine xipht(ia,ib,ixipht)
      implicit integer(a-z)
      save
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( srec=81, ssrec=62 )
      character*(srec) auxlin
      common/z22g/auxlin
      common/z30g/stib(1:sibuff)
      ii=0
      if(ib.gt.sibuff)then
      ii=1
      elseif(ia.le.0)then
      ii=1
      elseif(ib.lt.ia)then
      ii=1
      endif
      if(ixipht.lt.0)then
      if(ia+ixipht.le.0)then
      ii=1
      endif
      j1=ia
      j2=ib
      j3=1
      elseif(ixipht.gt.0)then
      if(ib+ixipht.gt.sibuff)then
      ii=1
      endif
      j1=ib
      j2=ia
      j3=-1
      endif
      if(ii.ne.0)then
      auxlin(1:srec)='xipht_1'
      call messag(1,0,0,0)
      endif
      if(ixipht.ne.0)then
      do i1=j1,j2,j3
      stib(i1+ixipht)=stib(i1)
      enddo
      endif
      return
      end
      subroutine cxipht(ia,ib,ixipht)
      implicit integer(a-z)
      save
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( srec=81, ssrec=62 )
      character*(srec) auxlin
      common/z22g/auxlin
      character*(scbuff) stcb
      common/z31g/stcb
      ii=0
      if(ia.le.0)then
      ii=1
      elseif(ib.lt.ia)then
      ii=1
      elseif(scbuff.lt.ib)then
      ii=1
      endif
      if(ixipht.lt.0)then
      if(ia+ixipht.le.0)then
      ii=1
      endif
      j1=ia
      j2=ib
      j3=1
      elseif(ixipht.gt.0)then
      if(ib+ixipht.gt.scbuff)then
      ii=1
      endif
      j1=ib
      j2=ia
      j3=-1
      endif
      if(ii.ne.0)then
      auxlin(1:srec)='cxipht_1'
      call messag(1,0,0,0)
      endif
      if(ixipht.ne.0)then
      do i1=j1,j2,j3
      stcb(i1+ixipht:i1+ixipht)=stcb(i1:i1)
      enddo
      endif
      return
      end
      subroutine trm(ia,ib)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      character*(srec) auxlin
      common/z22g/auxlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      if((ia.le.0).or.(ib.le.0))then
      auxlin(1:srec)='trm_1'
      call messag(1,0,0,0)
      endif
      ab=ia*ib
      st0=stibs(1)-ab+1
      do i1=1,ab-2
      ix=0
      jj=i1
   10 continue
      ii=jj/ia
      kk=ii+ib*(jj-ii*ia)
      if(ix.eq.1)then
      ytmp=stib(st0+kk)
      stib(st0+kk)=xtmp
      xtmp=ytmp
      endif
      if(kk.gt.i1)then
      jj=kk
      goto 10
      elseif(kk.eq.i1)then
      if((ix.eq.0).and.(kk.ne.jj))then
      ix=1
      jj=i1
      xtmp=stib(st0+i1)
      goto 10
      endif
      endif
      enddo
      return
      end
      subroutine qread(iunit,junit,nlin,slin)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( srecx=86, ich1=32, ich2=126 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( sxbuff=2040 )
      character*(scbuff) stcb
      common/z31g/stcb
      character*(srec) auxlin
      common/z22g/auxlin
      character*(sxbuff) stxb
      common/z27g/stxb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z52g/acf0(0:127),acf1(0:127),lfeed,aspace,dquote,squote
      ios=1
      read(unit=iunit,fmt=80,iostat=ios)stxb(1:srecx)
   80 format(a)
      nlin=nlin+1
      if(ios.ne.0)then
      if(ios.gt.0)then
      auxlin(1:srec)='run-time error while reading'
      call messag(1,0,0,junit)
      endif
      ios=1
      close(unit=iunit,status='keep',iostat=ios)
      if(ios.ne.0)then
      auxlin(1:srec)='system/filesystem error'
      call messag(1,0,0,0)
      endif
      iunit=0
      slin=-1
      goto 90
      endif
      do i1=srecx,1,-1
      if(ichar(stxb(i1:i1)).ne.aspace)then
      if(i1.lt.srec)then
      slin=i1
      goto 20
      else
      auxlin(1:srec)='line too long,'
      call messag(1,nlin,nlin,junit)
      endif
      endif
      enddo
      slin=0
      goto 40
   20 continue
      j1=0
      if(nlin.eq.1)then
      if(ichar(stxb(1:1)).eq.239)then
      if(ichar(stxb(2:2)).eq.187)then
      if(ichar(stxb(3:3)).eq.191)then
      j1=3
      slin=slin-j1
      endif
      endif
      endif
      endif
      do i1=j1+1,j1+slin
      ii=ichar(stxb(i1:i1))
      if((ii.lt.ich1).or.(ii.gt.ich2))then
      auxlin(1:srec)='wrong syntax,'
      call messag(1,nlin,nlin,junit)
      endif
      enddo
   40 continue
      ii=slin+1
      call vaocb(ii)
      jj=stcbs(1)
      do i1=jj+1,jj+slin
      j1=j1+1
      stcb(i1:i1)=stxb(j1:j1)
      enddo
      i1=jj+ii
      stcb(i1:i1)=char(lfeed)
   90 continue
      return
      end
      subroutine iopen(j1,j2,junit,jstat)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( lun1=7, lun2=63 )
      character*(srec) auxlin
      common/z22g/auxlin
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      logical lunit
      ii=0
      if(j1.le.0)then
      ii=1
      elseif(j2.le.0)then
      ii=1
      elseif(junit.le.0)then
      ii=1
      elseif(jstat.lt.0)then
      ii=1
      elseif(jstat.gt.1)then
      ii=1
      endif
      if(ii.ne.0)then
      auxlin(1:srec)='prefop_1'
      call messag(1,0,0,0)
      endif
      ios=1
      if(jstat.eq.0)then
      lunit=.false.
      else
      lunit=.true.
      endif
      inquire(file=stcb(j1:j1-1+j2),exist=lunit,iostat=ios)
      if(ios.ne.0)then
      auxlin(1:srec)='system/filesystem error'
      call messag(1,0,0,0)
      elseif(jstat.eq.0)then
      if(.not.lunit)then
      auxlin(1:srec)='could not be found'
      call messag(1,0,0,-junit)
      endif
      elseif(jstat.eq.1)then
      if(lunit)then
      auxlin(1:srec)='already exists'
      call messag(1,0,0,-junit)
      endif
      endif
      iunit=lun1
   11 continue
      if(iunit.gt.lun2)then
      auxlin(1:srec)='no logical unit number available'
      call messag(1,0,0,0)
      endif
      ios=1
      lunit=.true.
      inquire(unit=iunit,opened=lunit,iostat=ios)
      if(ios.ne.0)then
      auxlin(1:srec)='system/filesystem error'
      call messag(1,0,0,0)
      elseif(lunit)then
      iunit=iunit+1
      goto 11
      endif
      ios=1
      if(jstat.eq.0)then
      auxlin(1:3)='old'
      else
      auxlin(1:3)='new'
      endif
      open(unit=iunit,file=stcb(j1:j1-1+j2),status=auxlin(1:3),
     :iostat=ios,access='sequential')
      if(ios.ne.0)then
      auxlin(1:srec)='could not be opened'
      call messag(1,0,0,-junit)
      endif
      junit=iunit
      return
      end
      subroutine inputs
      implicit integer(a-z)
      save
      parameter ( maxleg=11, maxrho=7, maxi=11, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( srec=81, ssrec=62 )
      parameter ( eoa=-127, nap=-2047 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( sxbuff=2040 )
      common/z1in/leg(1:maxleg),nleg,nloop,incom
      common/z3in/momep(0:maxleg),momel(0:maxleg)
      common/z5in/dunit,munit,ounit,sunit,funit
      common/z6in/lffile,ffilea,ffileb
      common/z7in/lsfile,sfilea,sfileb
      common/z8in/lmfile,mfilea,mfileb
      common/z9in/lofile,ofilea,ofileb
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1
      common/z3g/nv(0:0),dpntro(0:0),vparto(0:0),vdeg(0:0),nvert
      common/z11g/nphi,nblok,nprop,npprop
      common/z12g/jflag(1:11),mflag(1:24)
      common/z14g/zcho(0:maxli),zbri(0:maxli),zpro(0:maxli),
     :rbri(0:maxli),sbri(0:maxli)
      common/z20g/tftyp(0:0),tfnarg(0:0),tfa(0:0),tfb(0:0),tfc(0:0),
     :tfo(0:0),tf2(0:0),ntf
      character*(srec) auxlin
      common/z22g/auxlin
      character*(sxbuff) stxb
      common/z27g/stxb
      common/z29g/pkey(0:0),wkey(0:0),prevl(0:0),nextl(0:0),popt3(0:0),
     :wopt3(0:0),fopt3(0:0),vopt3(0:0),popt5(0:0),wopt5(0:0),copt5(0:0)
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z33g/namep(0:0),namel(0:0)
      common/z34g/link(0:0),antiq(0:0),blok(0:0)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z37g/drecp(0:0),drecl(0:0),drecii(0:0),irecc(0:0),
     :frecc(0:0),ndrec,ncom
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      common/z42g/pmkr(0:0),pmkvma(0:0),pmkvmi(0:0)
      common/z43g/vmkr(0:0),vmkmao(0:0),vmkmio(0:0)
      common/z44g/xtstrp(0:0),xtstrl(0:0)
      common/z45g/qdatp,qdatl,qvp,qvl,dprefp,dprefl
      common/z46g/popt1(0:0),wopt1(0:0),copt1(0:0),popt7(0:0),
     :wopt7(0:0),copt7(0:0)
      common/z47g/ndiagp,ndiagl,hhp,hhl,doffp,doffl,noffp,noffl,wsint
      common/z49g/popt9(0:0),wopt9(0:0),copt9(0:0)
      common/z50g/tfta(0:0),tftb(0:0),tftic(0:0),ntft
      common/z51g/aplus,aminus,azero,anine,aeq
      common/z52g/acf0(0:127),acf1(0:127),lfeed,aspace,dquote,squote
      integer dreci(0:0),comii(0:0),comfi(0:0)
      integer xli(1:maxleg)
      datls=4
      incom=0
      outgo=0
      nleg=0
      dunit=1
      call iopen(qdatp,qdatl,dunit,0)
      ndrec=0
      newc=1
      nlin=0
      dbl=0
      level=0
      phase=1
   70 continue
      call qread(dunit,1,nlin,slin)
      if(slin.eq.-1)then
      if(newc.eq.0)then
      goto 521
      endif
      if(level.ne.-1)then
      if(level.lt.8)then
      auxlin(1:srec)='is incomplete'
      call messag(1,0,0,-1)
      endif
      endif
      goto 404
      endif
      jj=0
      if(slin.eq.0)then
      jj=1
      else
      ii=stcbs(1)+1
      if(stcb(ii:ii).eq.'#')then
      jj=1
      elseif(stcb(ii:ii).eq.'%')then
      jj=1
      elseif(stcb(ii:ii).eq.'*')then
      jj=1
      endif
      endif
      if(jj.ne.0)then
      if(newc.eq.1)then
      goto 70
      else
      goto 521
      endif
      endif
      if(newc.eq.1)then
      dbl=nlin
      bplic=1
      endif
      top=stcbs(1)+slin
      print *,stcb(stcbs(1)+1:top)
      jflag(7)=0
      do i1=top,stcbs(1)+1,-1
      j1=ichar(stcb(i1:i1))
      if(j1.ne.aspace)then
      if(j1.eq.squote)then
      bplic=-bplic
      endif
      bot=i1
      endif
      enddo
      if(bplic.ne.1)then
      goto 521
      endif
      call aoib(datls)
      ndrec=ndrec+1
      stib(stibs(1)-3)=stcbs(1)+1
      stib(stibs(1)-2)=slin+1
      stib(stibs(1)-1)=nlin
      stib(stibs(1))=nlin-dbl
      ii=slin+1
      call aocb(ii)
      if(stcb(top:top).eq.';')then
      newc=1
      else
      newc=0
      goto 70
      endif
      jj=stibs(1)-datls+1
      bot=stib(jj-datls*stib(stibs(1)))
      call stpa(stcb,bot,top,0)
      if(stib(stibs(1)+1).lt.3)then
      goto 521
      elseif(stib(stibs(1)+2).ne.3)then
      goto 521
      elseif(stib(stibs(1)+5).ne.1.or.stib(stibs(1)+6).ne.1)then
      goto 521
      endif
      j=stib(stibs(1)+3)
      jj=stib(stibs(1)+4)
   16 continue
      if(level.ge.0)then
      call mstr0(stcb,j,jj,pkey(0),wkey(0),ij)
      if(ij.eq.0)then
      if(level.gt.8)then
      level=-1
      goto 16
      else
      goto 521
      endif
      elseif(level.ne.stib(prevl(0)+ij))then
      goto 521
      endif
      if(level.gt.8)then
      jflag(4)=1
      endif
      level=stib(nextl(0)+ij)
      else
      call mstr0(stcb,j,jj,popt1(0),wopt1(0),ij)
      if(ij.eq.0)then
      goto 521
      endif
      if(stib(stibs(1)+8).ne.3)then
      goto 521
      endif
      i1=stib(stibs(1)+9)
      i2=stib(stibs(1)+10)
      call mstr0(stcb,i1,i2,popt5(0),wopt5(0),ij)
      if(ij.eq.0)then
      goto 521
      elseif(stib(copt5(0)+ij).lt.0)then
      goto 521
      endif
      endif
      goto 70
  404 continue
      ii=stibs(1)
      call aoib(datls)
      do i1=ii+1,stibs(1)
      stib(i1)=eoa
      enddo
      call trm(datls,ndrec+1)
      ii=ndrec+1
      drecii(0)=stibs(1)-ii
      dreci(0)=drecii(0)-ii
      drecl(0)=dreci(0)-ii
      drecp(0)=drecl(0)-ii
      ncom=0
      jcom=0
      j1=0
  556 continue
      j1=j1+1
      ii=stib(drecp(0)+j1)
      jj=ii-1
      j2=j1
  558 continue
      jj=jj+stib(drecl(0)+j2)
      if(j2.lt.ndrec)then
      if(stib(drecii(0)+j2+1).ne.0)then
      j2=j2+1
      goto 558
      endif
      endif
      j1=j2
      ncom=ncom+1
      call aoib(2)
      stib(stibs(1)-1)=ii
      stib(stibs(1))=jj
      if(j1.lt.ndrec)then
      goto 556
      endif
      call aoib(2)
      stib(stibs(1)-1)=eoa
      stib(stibs(1))=eoa
      call trm(2,ncom+1)
      ii=ncom+1
      comfi(0)=stibs(1)-ii
      comii(0)=comfi(0)-ii
      irecc(0)=stibs(1)
      call aoib(2*ii)
      frecc(0)=irecc(0)+ii
      stib(frecc(0))=eoa
      stib(stibs(1))=eoa
      j1=0
      j2=0
  137 continue
      if(j1.lt.ndrec)then
      j1=j1+1
      else
      goto 143
      endif
      if(stib(drecii(0)+j1).eq.0)then
      if(j2.gt.0)then
      stib(frecc(0)+j2)=j1-1
      endif
      j2=j2+1
      stib(irecc(0)+j2)=j1
      endif
      if(j2.le.ncom)then
      goto 137
      endif
  143 continue
      stib(frecc(0)+ncom)=ndrec
      if((j1.ne.ndrec).or.(j2.ne.ncom))then
      auxlin(1:srec)='inputs_1'
      call messag(1,0,0,0)
      endif
      phase=2
      icom=1
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),0)
      j=stibs(1)+8
      lofile=0
      ii=stib(stibs(1)+1)-3
      do 257 i=1,ii
      if(stib(j).ne.2)then
      goto 521
      endif
      if(i.eq.1)then
      ofilea=stib(j+1)
      endif
      if(i.eq.ii)then
      ofileb=stib(j+2)
      lofile=ofileb-ofilea+1
      endif
      j=j+3
  257 continue
      if(lofile.gt.0)then
      k=stds(stcb,ofilea,ofileb,2,1)
      if(k.lt.0)then
      goto 521
      elseif(k.eq.0)then
      lofile=0
      elseif(k.lt.lofile-2)then
      lofile=k
      ofilea=stcbs(1)+1
      call aocb(k)
      ofileb=stcbs(1)
      else
      lofile=lofile-2
      ofilea=ofilea+1
      ofileb=ofileb-1
      endif
      endif
      if(lofile.gt.0)then
      jflag(1)=1
      endif
      icom=2
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),0)
      j=stibs(1)+8
      lsfile=0
      ii=stib(stibs(1)+1)-3
      do 259 i=1,ii
      if(stib(j).ne.2)then
      goto 521
      endif
      if(i.eq.1)then
      sfilea=stib(j+1)
      endif
      if(i.eq.ii)then
      sfileb=stib(j+2)
      lsfile=sfileb-sfilea+1
      endif
      j=j+3
  259 continue
      if(lsfile.gt.0)then
      k=stds(stcb,sfilea,sfileb,2,1)
      if(k.lt.0)then
      goto 521
      elseif(k.eq.0)then
      lsfile=0
      elseif(k.lt.lsfile-2)then
      lsfile=k
      sfilea=stcbs(1)+1
      call aocb(k)
      sfileb=stcbs(1)
      else
      lsfile=lsfile-2
      sfilea=sfilea+1
      sfileb=sfileb-1
      endif
      endif
      if(jflag(1).ne.0)then
      if(lsfile.eq.0)then
      auxlin(1:srec)='name of style-file is missing,'
      call messag(1,0,0,1)
      endif
      call style
      endif
      lffile=0
      icom=3
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),0)
      if(stib(stibs(1)+1).lt.4)then
      goto 521
      endif
      j=stibs(1)+8
      lmfile=0
      ii=stib(stibs(1)+1)-3
      do 261 i=1,ii
      if(stib(j).ne.2)then
      goto 521
      endif
      if(i.eq.1)then
      mfilea=stib(j+1)
      endif
      if(i.eq.ii)then
      mfileb=stib(j+2)
      lmfile=mfileb-mfilea+1
      endif
      j=j+3
  261 continue
      if(lmfile.gt.0)then
      k=stds(stcb,mfilea,mfileb,2,1)
      if(k.lt.0)then
      goto 521
      elseif(k.eq.0)then
      lmfile=0
      elseif(k.lt.lmfile-2)then
      lmfile=k
      mfilea=stcbs(1)+1
      call aocb(k)
      mfileb=stcbs(1)
      else
      lmfile=lmfile-2
      mfilea=mfilea+1
      mfileb=mfileb-1
      endif
      endif
      if(lmfile.eq.0)then
      goto 521
      endif
      call model
      icom=4
  263 continue
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),0)
      if(stib(stibs(1)+1).eq.3)then
      nphia=0
      else
      if(stib(stibs(1)+11).eq.1.and.stib(stibs(1)+12).eq.5)then
      j=stib(stibs(1)+1)-2
      nphia=j/5
      if((nphia.lt.0).or.(j.ne.5*nphia))then
      goto 521
      endif
      if(icom.eq.4)then
      jflag(5)=1
      elseif(incom.eq.0)then
      jflag(5)=1
      elseif(jflag(5).eq.0)then
      goto 521
      endif
      j=stibs(1)+8
      do 267 i=1,nphia
      if(stib(j).ne.3)then
      goto 521
      elseif(stib(j+3).ne.1)then
      goto 521
      elseif(stib(j+4).ne.5)then
      goto 521
      elseif(stib(j+6).ne.3)then
      goto 521
      elseif(stib(j+9).ne.1)then
      goto 521
      elseif(stib(j+10).ne.6)then
      goto 521
      elseif(stib(j+12).ne.1)then
      goto 521
      elseif(stib(j+13).ne.2)then
      if(i.ne.nphia)then
      goto 521
      endif
      endif
      nleg=nleg+1
      if(nleg.gt.maxleg)then
      auxlin(1:srec)='too many legs'
      call messag(1,0,0,0)
      endif
      momep(nleg)=stib(j+7)
      momel(nleg)=stib(j+8)-momep(nleg)+1
      j=j+15
  267 continue
      else
      j=stib(stibs(1)+1)-2
      nphia=j/2
      if((nphia.lt.0).or.(j.ne.2*nphia))then
      goto 521
      endif
      if(icom.eq.4)then
      jflag(5)=0
      elseif(jflag(5).ne.0)then
      goto 521
      endif
      j=stibs(1)+8
      do 265 i=1,nphia
      if(stib(j).ne.3)then
      goto 521
      elseif(stib(j+3).ne.1)then
      goto 521
      elseif(stib(j+4).ne.2)then
      if(i.ne.nphia)then
      goto 521
      endif
      endif
      nleg=nleg+1
      if(nleg.gt.maxleg)then
      auxlin(1:srec)='too many legs'
      call messag(1,0,0,0)
      endif
      call aocb(1)
      k=stcbs(1)
      if(icom.eq.4)then
      stcb(k:k)=char(112)
      else
      stcb(k:k)=char(113)
      endif
      momep(nleg)=k
      call dkar(i,k,jk)
      k=k+jk
      call aocb(k-stcbs(1))
      momel(nleg)=stcbs(1)-momep(nleg)+1
      j=j+6
  265 continue
      endif
      endif
      if(icom.eq.4)then
      ij=0
      else
      ij=incom
      endif
      ii=6+9*jflag(5)
      jj=stibs(1)+8-ii
      do 273 i=1,nphia
      jj=jj+ii
      do 271 i1=1,nphi
      ik=stib(namel(0)+i1)
      j=stib(jj+1)
      k=stib(jj+2)
      if(ik-1.eq.k-j)then
      jk=stib(namep(0)+i1)
      if(stcb(j:k).eq.stcb(jk:jk-1+ik))then
      ij=ij+1
      if(icom.eq.4)then
      leg(ij)=stib(link(0)+i1)
      else
      leg(ij)=i1
      endif
      goto 273
      endif
      endif
  271 continue
      auxlin(1:srec)='unknown external particle(s)'
      call messag(1,0,0,0)
  273 continue
      if(icom.eq.4)then
      icom=5
      incom=nphia
      goto 263
      else
      outgo=nphia
      endif
      if(nleg.ne.incom+outgo)then
      auxlin(1:srec)='inputs_2'
      call messag(1,0,0,0)
      endif
      if(nleg.eq.0)then
      auxlin(1:srec)='no external particles listed'
      call messag(1,0,0,0)
      endif
      icom=6
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),0)
      if(stib(stibs(1)+1).ne.4)then
      goto 521
      elseif(stib(stibs(1)+8).ne.4)then
      goto 521
      endif
      nloop=stoz(stcb,stib(stibs(1)+9),stib(stibs(1)+10))
      if(nloop.lt.0)then
      auxlin(1:srec)='check loops in'
      call messag(1,0,0,1)
      elseif(nloop.gt.maxrho.or.nloop+nleg.gt.max(maxleg,maxrho))then
      auxlin(1:srec)='too many legs and/or loops'
      call messag(1,0,0,0)
      endif
      if(nleg+nloop.lt.2)then
      auxlin(1:srec)='check legs, loops in'
      call messag(1,0,0,1)
      endif
      if(nleg.eq.2.and.nloop.eq.0)then
      auxlin(1:srec)='case legs=2, loops=0 not accepted'
      call messag(1,0,0,0)
      endif
      icom=7
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),0)
      if(stib(stibs(1)+1).eq.3)then
      call aocb(1)
      stcb(stcbs(1):stcbs(1))=char(107)
      momep(0)=stcbs(1)
      momel(0)=1
      elseif(stib(stibs(1)+1).eq.4)then
      if(stib(stibs(1)+8).ne.3)then
      goto 521
      endif
      momep(0)=stib(stibs(1)+9)
      momel(0)=stib(stibs(1)+10)-stib(stibs(1)+9)+1
      else
      goto 521
      endif
      if(jflag(1).ne.0)then
      call vaocb(momel(0))
      k=stcbs(1)+momel(0)
      stcb(stcbs(1)+1:k)=stcb(momep(0):momep(0)-1+momel(0))
      do 289 i=1,nloop
      ik=k
      call dkar(i,ik,jk)
      ik=ik+jk
      do 287 j=1,nleg
      if(ik-momel(j).eq.k-momel(0))then
      if(stcb(stcbs(1)+1:ik).eq.
     :stcb(momep(j):momep(j)-1+momel(j)))then
      auxlin(1:srec)='conflict between names of external and '//
     :'internal momenta'
      call messag(1,0,0,0)
      endif
      endif
  287 continue
  289 continue
      endif
      icom=8
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),0)
      jj=stib(stibs(1)+1)-2
      if(jj.eq.1)then
      ii=0
      else
      ii=jj/2
      if(jj.ne.2*ii)then
      goto 521
      endif
      endif
      do i1=1,ii
      jj=stibs(1)+6*i1
      if((stib(jj+5).ne.1).or.(stib(jj+6).ne.2))then
      if(i1.ne.ii)then
      goto 521
      endif
      endif
      j1=stib(jj+3)
      j2=stib(jj+4)
      call mstr0(stcb,j1,j2,popt3(0),wopt3(0),ij)
      if(ij.ne.0)then
      j1=stib(fopt3(0)+ij)
      j2=stib(vopt3(0)+ij)
      if(mflag(j1).ne.0)then
      call wput(11,0,0)
      if(mflag(j1)*j2.lt.0)then
      jflag(2)=1
      elseif(mflag(j1).lt.j2)then
      mflag(j1)=j2
      endif
      endif
      mflag(j1)=j2
      endif
      if((ij.eq.0).or.(j1.eq.7))then
      auxlin(1:srec)='unknown option,'
      goto 522
      endif
      enddo
      if(jflag(4).ne.0)then
      icom=9
      call stpa(stcb,stib(comii(0)+icom),stib(comfi(0)+icom),0)
      jj=stibs(1)+1
      if(stib(jj).ne.4)then
      goto 521
      endif
      jj=stibs(1)+8
      if(stib(jj).ne.4)then
      goto 521
      endif
      jj=jj+1
      j1=stib(jj)
      j2=stib(jj+1)
      if(j2.lt.j1)then
      goto 521
      elseif(stdz(stcb,j1,j2).eq.0)then
      goto 521
      endif
      if(nullz(stcb,j1,j2).eq.0)then
      if(ichar(stcb(j1:j1)).eq.aminus)then
      goto 521
      endif
      if(ichar(stcb(j1:j1)).eq.aplus)then
      j1=j1+1
      endif
      j3=j1-1
      ii=0
      do i1=j1,j2-1
      if(ii.eq.0)then
      if(ichar(stcb(i1:i1)).eq.azero)then
      j3=i1
      else
      ii=1
      endif
      endif
      enddo
      noffl=j2-j3
      if(noffl.ge.wsint)then
      ii=stoz(stcb,j1,j2)
      endif
      stcb(noffp+1:noffp+noffl)=stcb(j3+1:j2)
      doffl=noffl
      stcb(doffp+1:doffp+doffl)=stcb(j3+1:j2)
      endif
      endif
      ipass=1
      nphib=0
      ntf1=0
      ntf=0
      do i1=0,maxli
      zpro(i1)=1
      zbri(i1)=1
      rbri(i1)=1
      sbri(i1)=1
      if(i1.lt.nloop)then
      zcho(i1)=0
      else
      zcho(i1)=1
      endif
      enddo
   54 continue
      do 71 jcom=icom+1,ncom
      call stpa(stcb,stib(comii(0)+jcom),stib(comfi(0)+jcom),0)
      ii=stib(stibs(1)+1)
      ik=stibs(1)+3*ii
      if(ii.lt.7)then
      goto 521
      elseif(mod(ii,2).eq.0)then
      goto 521
      elseif(stib(stibs(1)+8).ne.3)then
      goto 521
      elseif(stib(stibs(1)+11).ne.1)then
      goto 521
      elseif(stib(stibs(1)+12).ne.5)then
      goto 521
      elseif(stib(ik-7).ne.4)then
      goto 521
      elseif(stib(ik-4).ne.1)then
      goto 521
      elseif(stib(ik-3).ne.6)then
      goto 521
      endif
      j1=stib(stibs(1)+3)
      j2=stib(stibs(1)+4)
      call mstr0(stcb,j1,j2,popt1(0),wopt1(0),ij)
      if(ij.le.0)then
      goto 521
      endif
      tfv=stib(copt1(0)+ij)
      j1=stib(stibs(1)+9)
      j2=stib(stibs(1)+10)
      call mstr0(stcb,j1,j2,popt5(0),wopt5(0),ij)
      if(ij.le.0)then
      goto 521
      endif
      atyp=stib(copt5(0)+ij)
      styp=atyp*tfv
      if(atyp.lt.12)then
      if(ii.lt.9)then
      goto 521
      elseif(stib(ik-13).ne.4)then
      goto 521
      elseif(stib(ik-10).ne.1)then
      goto 521
      elseif(stib(ik-9).ne.2)then
      goto 521
      endif
      endif
      if(atyp.eq.6)then
      mflag(20)=1
      elseif(atyp.eq.7)then
      mflag(21)=1
      elseif(atyp.eq.11)then
      mflag(22)=1
      elseif(atyp.eq.12)then
      mflag(23)=1
      elseif(atyp.gt.12)then
      goto 521
      endif
      j1=stib(stibs(1)+1)
      if(atyp.lt.11)then
      j3=max(0,(j1-9)/2)
      elseif(atyp.eq.11)then
      j3=max(0,(j1-11)/2)
      elseif(atyp.eq.12)then
      j3=max(0,(j1-5)/2)
      endif
      if(atyp.ge.11)then
      if(j3.le.0)then
      goto 521
      endif
      endif
      if(ipass.eq.1)then
      if(j3.ne.0)then
      ntf1=ntf1+1
      nphib=nphib+j3
      endif
      goto 71
      endif
      if(atyp.lt.12)then
      j1=stoz(stcb,stib(ik-12),stib(ik-11))
      j2=stoz(stcb,stib(ik-6),stib(ik-5))
      if(j1.gt.j2)then
      goto 521
      endif
      if(atyp.lt.6)then
      if(j1.lt.0)then
      goto 521
      endif
      if(atyp.gt.1)then
      jflag(10)=1
      endif
      elseif(atyp.eq.11)then
      if(j1.le.0)then
      goto 521
      endif
      endif
      else
      j1=0
      j2=0
      endif
      if(j3.eq.0)then
      goto 73
      endif
      ntf=ntf+1
      stib(tftyp(0)+ntf)=styp
      stib(tfa(0)+ntf)=j1
      stib(tfb(0)+ntf)=j2
      stib(tfc(0)+ntf)=0
      stib(tfnarg(0)+ntf)=j3
      stib(stib(tfo(0)+ntf))=eoa
      if(ntf.lt.ntf1)then
      stib(tfo(0)+ntf+1)=stib(tfo(0)+ntf)+stib(tfnarg(0)+ntf)+1
      endif
      if(atyp.gt.5)then
      goto 37
      endif
      do i1=1,stib(tfnarg(0)+ntf)
      j=stibs(1)+6*i1
      if(stib(j+8).ne.3)then
      goto 521
      elseif(stib(j+11).ne.1)then
      goto 521
      elseif(stib(j+12).ne.2)then
      goto 521
      endif
      j1=stib(j+9)
      j2=stib(j+10)
      call mstr0(stcb,j1,j2,namep(0),namel(0),jk)
      if(jk.eq.0)then
      auxlin(1:srec)='unknown field,'
      goto 522
      endif
      stib(stib(tfo(0)+ntf)+i1)=jk
      enddo
      i4=stib(tfo(0)+ntf)
      call hsort(stib,i4+1,i4+stib(tfnarg(0)+ntf))
      tfcl=stib(tfnarg(0)+ntf)
      tfch=1
      tfcn=1
      i3=1
      do i1=2,stib(tfnarg(0)+ntf)
      if(stib(i4+i1).eq.stib(i4+i1-1).or.
     :stib(i4+i1).eq.stib(link(0)+stib(i4+i1-1)))then
      i3=i3+1
      if(i3.gt.tfch)then
      tfch=i3
      endif
      else
      if(i3.lt.tfcl)then
      tfcl=i3
      endif
      tfcn=tfcn+1
      i3=1
      endif
      enddo
      if(i3.lt.tfcl)then
      tfcl=i3
      endif
      if(tfcn.ne.nprop)then
      tfcl=0
      endif
      if(styp.gt.0)then
      do i2=0,maxli
      if(i2*tfch.lt.stib(tfa(0)+ntf).or.
     :(i2*tfcl.gt.stib(tfb(0)+ntf)))then
      if(styp.eq.1)then
      zpro(i2)=0
      elseif(styp.eq.2)then
      zbri(i2)=0
      elseif(styp.eq.3)then
      zcho(i2)=0
      elseif(styp.eq.4)then
      rbri(i2)=0
      elseif(styp.eq.5)then
      sbri(i2)=0
      endif
      endif
      enddo
      elseif(styp.lt.0)then
      do i2=0,maxli
      if(i2*tfch.le.stib(tfb(0)+ntf))then
      if(i2*tfcl.ge.stib(tfa(0)+ntf))then
      if(styp.eq.-1)then
      zpro(i2)=0
      elseif(styp.eq.-2)then
      zbri(i2)=0
      elseif(styp.eq.-3)then
      zcho(i2)=0
      elseif(styp.eq.-4)then
      rbri(i2)=0
      elseif(styp.eq.-5)then
      sbri(i2)=0
      endif
      endif
      endif
      enddo
      endif
      goto 71
   73 continue
      if(styp.gt.0)then
      do i2=0,maxli
      if((i2.lt.j1).or.(i2.gt.j2))then
      if(styp.eq.1)then
      zpro(i2)=0
      elseif(styp.eq.2)then
      zbri(i2)=0
      elseif(styp.eq.3)then
      zcho(i2)=0
      elseif(styp.eq.4)then
      rbri(i2)=0
      elseif(styp.eq.5)then
      sbri(i2)=0
      endif
      endif
      enddo
      else
      do i2=0,maxli
      if(i2.ge.j1)then
      if(i2.le.j2)then
      if(styp.eq.-1)then
      zpro(i2)=0
      elseif(styp.eq.-2)then
      zbri(i2)=0
      elseif(styp.eq.-3)then
      zcho(i2)=0
      elseif(styp.eq.-4)then
      rbri(i2)=0
      elseif(styp.eq.-5)then
      sbri(i2)=0
      endif
      endif
      endif
      enddo
      endif
      goto 71
   37 continue
      if(atyp.eq.6)then
      if(stib(tfnarg(0)+ntf).ne.1)then
      goto 521
      endif
      j1=stib(stibs(1)+15)
      j2=stib(stibs(1)+16)
      call mstr0(stcb,j1,j2,vmkp(0),vmkl(0),kid)
      if(kid.eq.0)then
      auxlin(1:srec)='wrong argument in vsum,'
      goto 522
      endif
      if(stib(vmkr(0)+kid).ne.4)then
      auxlin(1:srec)='invalid keyword in vsum,'
      goto 522
      endif
      stib(vmks(0)+kid)=1
      stib(stib(tfo(0)+ntf)+1)=kid
      elseif(atyp.eq.7)then
      if(stib(tfnarg(0)+ntf).ne.1)then
      goto 521
      endif
      j1=stib(stibs(1)+15)
      j2=stib(stibs(1)+16)
      call mstr0(stcb,j1,j2,pmkp(0),pmkl(0),kid)
      if(kid.eq.0)then
      auxlin(1:srec)='wrong argument in psum,'
      goto 522
      endif
      if((stib(pmkr(0)+kid).ne.4).or.(stib(pmkd(0)+kid).ne.1))then
      auxlin(1:srec)='invalid keyword in psum,'
      goto 522
      endif
      stib(stib(tfo(0)+ntf)+1)=kid
      j1=stib(pmkvmi(0)+kid)
      j2=stib(pmkvma(0)+kid)
      do i2=0,maxli
      if(styp.gt.0)then
      if(i2*j1.gt.stib(tfb(0)+ntf).or.
     :i2*j2.lt.stib(tfa(0)+ntf))then
      zpro(i2)=0
      endif
      else
      if(i2*j1.ge.stib(tfa(0)+ntf).and.
     :i2*j2.le.stib(tfb(0)+ntf))then
      zpro(i2)=0
      endif
      endif
      enddo
      elseif((atyp.eq.11).or.(atyp.eq.12))then
      j3=stib(tfnarg(0)+ntf)
      if(atyp.eq.11)then
      if(j3.lt.stib(tfb(0)+ntf))then
      goto 521
      endif
      elseif(atyp.eq.12)then
      if(nleg.lt.2)then
      call wput(-14,0,0)
      endif
      if(j3.ge.nleg)then
      goto 521
      endif
      endif
      do i1=1,nleg
      xli(i1)=0
      enddo
      do i1=1,j3
      j1=stibs(1)+9+6*i1
      if(stib(j1-1).ne.4)then
      goto 521
      endif
      j2=stoz(stcb,stib(j1),stib(j1+1))
      if(j2.ge.0)then
      goto 521
      elseif(j2+2*nleg.lt.0)then
      goto 521
      endif
      j1=(-j2)/2
      if(j2+2*j1.eq.0)then
      j1=j1+incom
      if(j1.gt.nleg)then
      goto 521
      endif
      else
      j1=j1+1
      if(j1.gt.incom)then
      goto 521
      endif
      endif
      if(xli(j1).ne.0)then
      goto 521
      endif
      xli(j1)=1
      stib(stib(tfo(0)+ntf)+i1)=j1
      enddo
      if(atyp.eq.11)then
      j1=stibs(1)+15+6*j3
      call mstr0(stcb,stib(j1),stib(j1+1),popt9(0),wopt9(0),ij)
      if(ij.eq.0)then
      auxlin(1:srec)='invalid elink statement,'
      goto 522
      endif
      j1=stib(copt9(0)+ij)
      stib(tfc(0)+ntf)=j1
      j2=0
      if(j1.ne.0)then
      if(j3.eq.nleg)then
      if(stib(tfa(0)+ntf).eq.1)then
      if(stib(tfb(0)+ntf).eq.nleg)then
      j2=styp
      endif
      endif
      endif
      else
      if(stib(tfa(0)+ntf).eq.1)then
      if(stib(tfb(0)+ntf).eq.j3)then
      j2=styp
      endif
      endif
      endif
      if(j2.ne.0)then
      if(j2.gt.0)then
      call wput(12,0,0)
      else
      jflag(2)=1
      call wput(13,0,0)
      endif
      endif
      endif
      else
      auxlin(1:srec)='inputs_3'
      call messag(1,0,0,0)
      endif
   71 continue
      if(ipass.eq.1)then
      ii=ntf1+1
      jj=7*ii+nphib+ntf1
      tftyp(0)=stibs(1)
      call aoib(jj)
      tfnarg(0)=tftyp(0)+ii
      tfa(0)=tfnarg(0)+ii
      tfb(0)=tfa(0)+ii
      tfc(0)=tfb(0)+ii
      tf2(0)=tfc(0)+ii
      tfo(0)=tf2(0)+ii
      stib(tfo(0)+1)=tfo(0)+ii
      stib(tfnarg(0))=eoa
      stib(tfa(0))=eoa
      stib(tfb(0))=eoa
      stib(tfc(0))=eoa
      stib(tf2(0))=eoa
      stib(tfo(0))=eoa
      stib(stibs(1))=eoa
      ipass=2
      goto 54
      endif
      call spp(1)
      call vsig
      do i1=1,ntft
      stib(tfta(0)+i1)=0
      stib(tftb(0)+i1)=0
      enddo
      j1=tftic(0)
      j3=1
      jj=0
   29 continue
      j1=j1+1
      j2=stib(j1)
      if(j2.ge.0)then
      if(j2.gt.0)then
      j3=j2
      jj=jj+1
      endif
      goto 29
      endif
      j1=j1-tftic(0)-1
      ii=0
      do i1=1,j1
      j2=stib(tftic(0)+i1)
      if(j2.gt.0)then
      do i2=1,ntf
      if(abs(stib(tftyp(0)+i2)).eq.i1)then
      ii=ii+1
      j3=stib(tfta(0)+j2)
      if(j3.eq.0)then
      stib(tfta(0)+j2)=ii
      endif
      stib(tftb(0)+j2)=ii
      stib(tf2(0)+ii)=i2
      endif
      enddo
      endif
      enddo
      do i1=1,ntft
      if(stib(tfta(0)+i1).eq.0)then
      stib(tfta(0)+i1)=1
      endif
      enddo
      ii=1
      do i1=1,nleg
      if(stib(antiq(0)+leg(i1)).ne.0)then
      ii=-ii
      endif
      enddo
      if(ii.lt.0)then
      call wput(10,0,0)
      jflag(2)=1
      endif
      if(nleg.gt.1)then
      ii=0
      jj=stib(blok(0)+stib(link(0)+leg(1)))
      do i1=2,nleg
      if(ii.eq.0)then
      if(stib(blok(0)+stib(link(0)+leg(i1))).ne.jj)then
      ii=1
      call wput(9,0,0)
      jflag(2)=1
      endif
      endif
      enddo
      endif
      if((nleg.ne.2).or.(nloop.ne.0))then
      do i1=1,nleg
      jj=stib(link(0)+leg(i1))
      do i2=1,nrho
      if(stib(nv(0)+i2).gt.0)then
      if(stib(stib(stib(dpntro(0)+i2)+jj)+1).eq.jj)then
      goto 400
      endif
      endif
      enddo
      jflag(2)=1
      kk=0
      if(i1.gt.incom)then
      ii=stib(namep(0)+leg(i1))
      do i2=incom+1,i1-1
      if(leg(i2).eq.leg(i1))then
      kk=1
      endif
      enddo
      if(kk.eq.0)then
      j1=ii-1+stib(namel(0)+leg(i1))
      call wput(8,ii,j1)
      endif
      else
      ii=stib(namep(0)+jj)
      do i2=1,i1-1
      if(leg(i2).eq.leg(i1))then
      kk=1
      endif
      enddo
      if(kk.eq.0)then
      j1=ii-1+stib(namel(0)+jj)
      call wput(7,ii,j1)
      endif
      endif
  400 continue
      enddo
      endif
      if(nphi.eq.1)then
      mflag(14)=1
      endif
      ii=0
      jj=0
      do i1=1,nphi
      if(stib(antiq(0)+i1).eq.0)then
      ii=1
      else
      jj=1
      endif
      enddo
      if(ii.eq.0)then
      mflag(15)=-1
      elseif(jj.eq.0)then
      mflag(15)=1
      endif
      mflag(16)=1
      j1=0
   02 continue
      if(j1.lt.nvert)then
      j1=j1+1
      j2=stib(vparto(0)+j1)
      ii=-1
      do i2=1,stib(vdeg(0)+j1)
      ii=ii+stib(antiq(0)+stib(j2+i2))
      enddo
      if(abs(ii).eq.1)then
      goto 02
      endif
      mflag(16)=0
      endif
      if(mflag(14).eq.0)then
      jflag(11)=0
      elseif(mflag(17).ne.0)then
      jflag(11)=0
      else
      jflag(11)=1
      endif
      goto 90
  521 continue
      auxlin(1:srec)='wrong syntax,'
  522 continue
      if(phase.eq.1)then
      call messag(1,dbl,nlin,1)
      else
      j1=1
      j2=0
      if(jcom.eq.0)then
      jcom=icom
      endif
  523 continue
      if(stib(drecii(0)+j1).eq.0)then
      j2=j2+1
      nlin=stib(dreci(0)+j1)
      endif
      if(j1.lt.ndrec.and.j2.lt.jcom)then
      j1=j1+1
      goto 523
      endif
  524 continue
      if(j1.lt.ndrec)then
      if(stib(drecii(0)+j1+1).ne.0)then
      j1=j1+1
      goto 524
      endif
      endif
      call messag(1,nlin,nlin+stib(drecii(0)+j1),1)
      endif
   90 continue
      return
      end
      subroutine model
      implicit integer(a-z)
      save
      parameter ( maxleg=11, maxrho=7, maxi=11, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( srec=81, ssrec=62 )
      parameter ( eoa=-127, nap=-2047 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      common/z1in/leg(1:maxleg),nleg,nloop,incom
      common/z5in/dunit,munit,ounit,sunit,funit
      common/z6in/lffile,ffilea,ffileb
      common/z8in/lmfile,mfilea,mfileb
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1
      common/z3g/nv(0:0),dpntro(0:0),vparto(0:0),vdeg(0:0),nvert
      common/z11g/nphi,nblok,nprop,npprop
      character*(srec) auxlin
      common/z22g/auxlin
      common/z29g/pkey(0:0),wkey(0:0),prevl(0:0),nextl(0:0),popt3(0:0),
     :wopt3(0:0),fopt3(0:0),vopt3(0:0),popt5(0:0),wopt5(0:0),copt5(0:0)
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z32g/sucpal(0:11,0:11)
      common/z33g/namep(0:0),namel(0:0)
      common/z34g/link(0:0),antiq(0:0),blok(0:0)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      common/z38g/gmkp(0:0),gmkl(0:0),gmkd(0:0),gmko(0:0),
     :gmkvp(0:0),gmkvl(0:0)
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      common/z46g/popt1(0:0),wopt1(0:0),copt1(0:0),popt7(0:0),
     :wopt7(0:0),copt7(0:0)
      common/z51g/aplus,aminus,azero,anine,aeq
      common/z52g/acf0(0:127),acf1(0:127),lfeed,aspace,dquote,squote
      integer apmkd(0:0),avmkd(0:0),llp(0:0)
      if(lffile.gt.0)then
      xfile=0
      else
      xfile=1
      endif
   05 continue
      xfile=xfile+1
      if(xfile.eq.1)then
      funit=3
      call iopen(ffilea,lffile,funit,0)
      elseif(xfile.eq.2)then
      munit=2
      call iopen(mfilea,lmfile,munit,0)
      endif
      nlin=0
      dbl=0
      newc=1
      level=1
      if(xfile.eq.1)then
      lastp=0
      lastkp=0
      else
      llp(0)=nap
      ngmk=0
      npmk=0
      nvmk=0
      nrho=0
      nphi=0
      nprop=0
      nvert=0
      endif
   65 continue
      dip=0
      nest=0
   70 continue
      if(xfile.eq.1)then
      goto 05
      elseif(xfile.eq.2)then
      call qread(munit,2,nlin,slin)
      if(slin.eq.-1)then
      goto 777
      endif
      endif
      if(newc.eq.1)then
      dbl=nlin
      endif
      jj=0
      ii=stcbs(1)+1
      if(slin.eq.0)then
      jj=1
      elseif(stcb(ii:ii).eq.'#')then
      jj=1
      elseif(stcb(ii:ii).eq.'%')then
      jj=1
      elseif(stcb(ii:ii).eq.'*')then
      jj=1
      endif
      if(jj.ne.0)then
      if(newc.eq.1)then
      goto 70
      else
      goto 521
      endif
      endif
      ibot=stcbs(1)+1
      itop=stcbs(1)+slin
      do i1=ibot,itop
      if(ichar(stcb(i1:i1)).ne.aspace)then
      ibot=i1
      goto 16
      endif
      enddo
   16 continue
      ii=0
      quote=0
      do i1=ibot,itop
      jj=ichar(stcb(i1:i1))
      if(jj.eq.squote)then
      if(nest.ne.1)then
      goto 521
      endif
      if(quote.eq.0)then
      quote=squote
      endif
      if(quote.eq.squote)then
      ii=1-ii
      endif
      elseif(jj.eq.dquote)then
      if(nest.ne.1)then
      goto 521
      endif
      if(quote.eq.0)then
      quote=dquote
      endif
      if(quote.eq.dquote)then
      ii=1-ii
      endif
      elseif(ii.eq.0)then
      if(stcb(i1:i1).eq.'[')then
      if(nest.ne.0)then
      goto 521
      endif
      nest=1
      elseif(stcb(i1:i1).eq.']')then
      if(nest.ne.1)then
      goto 521
      endif
      nest=2
      elseif(stcb(i1:i1).eq.',')then
      if(nest.ne.1)then
      goto 521
      endif
      quote=0
      elseif(jj.ne.aspace)then
      if(nest.ne.1)then
      goto 521
      endif
      endif
      else
      if(nest.ne.1)then
      goto 521
      endif
      endif
      enddo
      if(ii.ne.0)then
      goto 521
      endif
      if(dbl.eq.nlin)then
      jbot=ibot
      endif
      jj=slin+1
      call aocb(jj)
      if(stcb(itop:itop).ne.']')then
      newc=0
      goto 70
      endif
      if(nest.ne.2)then
      goto 521
      endif
      newc=1
      ibot=jbot
      top=stcbs(1)
   50 continue
      pal=0
      mpal=11
      vel=0
      vin=0
      idin=0
      kin=0
      kid=0
      knf=0
      id1=0
      id2=0
      is1=0
      is2=0
      plm=0
      ix=ibot
  280 continue
      ic=ichar(stcb(ix:ix))
      if(sucpal(pal,acf0(ic)).lt.0)then
      goto 521
      endif
      if(pal.lt.3)then
      if(pal.eq.1)then
      if(acf0(ic).eq.0)then
      if(id1.eq.0)then
      id1=ix
      endif
      id2=ix
      elseif(acf0(ic).gt.2)then
      if(stcb(ix:ix).eq.',')then
      if(level.eq.1)then
      if(ngmk.eq.0)then
      call aoib(1)
      stib(stibs(1))=eoa
      llp(0)=stibs(1)
      endif
      call rgmki(llp(0))
      if(xfile.eq.1)then
      goto 521
      endif
      level=2
      fpass=1
      call aoib(1)
      stib(stibs(1))=eoa
      llp(0)=stibs(1)
      lastkp=0
      off1=0
      endif
      elseif(stcb(ix:ix).eq.';')then
      if(level.lt.2)then
      goto 521
      endif
      elseif(stcb(ix:ix).eq.'=')then
      if(level.ne.1.or.kin.ne.0)then
      goto 521
      endif
      kin=1
      goto 283
      elseif(stcb(ix:ix).eq.']')then
      if(level.lt.2)then
      goto 521
      endif
      endif
      if(level.eq.2)then
      idin=idin+1
      if(idin.gt.4)then
      goto 521
      endif
      goto 281
      elseif(level.eq.3)then
      idin=idin+1
      goto 282
      endif
      endif
      elseif(pal.eq.2)then
      if(stcb(ix:ix).eq.'=')then
      kin=kin+1
      goto 283
      elseif(acf0(ic).eq.0)then
      if(id1.eq.0)then
      id1=ix
      endif
      id2=ix
      endif
      endif
      elseif(pal.lt.6)then
      if(pal.eq.3)then
      if(vel.ne.0.or.vin.ne.0)then
      goto 521
      endif
      if(ichar(stcb(ix:ix)).eq.squote)then
      plm=1-plm
      if(is1.eq.0)then
      is1=ix
      endif
      is2=ix
      else
      if(acf0(ic).eq.0)then
      if(is1.eq.0)then
      is1=ix
      endif
      is2=ix
      elseif(stcb(ix:ix).eq.'(')then
      if(is1.ne.0)then
      goto 521
      endif
      vel=1
      elseif(stcb(ix:ix).eq.','.or.stcb(ix:ix).eq.']')then
      goto 285
      endif
      endif
      elseif(pal.eq.4)then
      elseif(pal.eq.5)then
      if(vel.ne.0.or.vin.ne.0.or.plm.ne.1)then
      goto 521
      endif
      if(ichar(stcb(ix:ix)).eq.squote)then
      plm=1-plm
      endif
      is2=ix
      endif
      elseif(pal.lt.9)then
      if(pal.eq.6)then
      if(vel.ne.1.or.plm.ne.0)then
      goto 521
      endif
      if(ichar(stcb(ix:ix)).eq.squote)then
      plm=1-plm
      if(is1.eq.0)then
      is1=ix
      endif
      is2=ix
      else
      if(acf0(ic).eq.0)then
      if(is1.eq.0)then
      is1=ix
      endif
      is2=ix
      elseif(stcb(ix:ix).eq.','.or.stcb(ix:ix).eq.')')then
      vin=vin+1
      goto 285
      endif
      endif
      elseif(pal.eq.7)then
      elseif(pal.eq.8)then
      if(vel.ne.1.or.plm.ne.1)then
      goto 521
      endif
      if(ichar(stcb(ix:ix)).eq.squote)then
      plm=1-plm
      endif
      is2=ix
      endif
      elseif(pal.lt.12)then
      if(pal.eq.9)then
      if(level.gt.2.or.level.le.0)then
      goto 521
      endif
      if(stcb(ix:ix).eq.','.or.stcb(ix:ix).eq.']')then
      vin=0
      vel=0
      endif
      elseif(pal.eq.10)then
      elseif(pal.eq.11)then
      endif
      else
      goto 521
      endif
  287 continue
      pal=sucpal(pal,acf0(ic))
      if(ix.lt.top)then
      ix=ix+1
      goto 280
      endif
      if(vel.ne.0)then
      goto 521
      elseif(pal.ne.mpal)then
      goto 521
      elseif(plm.ne.0)then
      goto 521
      elseif(id1.ne.0)then
      goto 521
      elseif(is1.ne.0)then
      goto 521
      endif
      if(level.eq.1)then
      goto 65
      elseif(level.eq.2)then
      if(fpass.eq.1)then
      else
      do i1=1,npmk
      if(stib(apmkd(0)+i1).eq.3)then
      if(stib(lastp+1).eq.0)then
      auxlin(1:srec)='wrong dimension for M-function,'
      call messag(1,dbl,nlin,4-xfile)
      endif
      endif
      if(stib(apmkd(0)+i1).ne.stib(pmkd(0)+i1))then
      ii=1
      if(stib(apmkd(0)+i1).eq.2)then
      if(stib(pmkd(0)+i1).eq.3)then
      ii=stib(lastp+1)
      endif
      endif
      if(ii.ne.0)then
      if(stib(apmkd(0)+i1).eq.0)then
      auxlin(1:srec)='missing M-function,'
      else
      auxlin(1:srec)='wrong dimension for M-function,'
      endif
      call messag(1,dbl,nlin,2)
      endif
      endif
      enddo
      nprop=nprop+1
      endif
      elseif(level.eq.3)then
      if(fpass.eq.0)then
      do i1=1,nvmk
      if(stib(avmkd(0)+i1).eq.0)then
      auxlin(1:srec)='missing M-function,'
      call messag(1,dbl,nlin,2)
      endif
      enddo
      jj=stib(lastp+1)
      if(jj.lt.3)then
      auxlin(1:srec)='vertex degree is too small,'
      call messag(1,dbl,nlin,2)
      endif
      ii=0
      do i1=1,jj
      if(stib(antiq(0)+stib(lastp+1+i1)).ne.0)then
      ii=1-ii
      endif
      enddo
      if(ii.ne.0)then
      auxlin(1:srec)='odd vertex,'
      call messag(1,dbl,nlin,2)
      endif
      endif
      endif
      if(level.eq.2)then
      if(fpass.eq.1)then
      ii=stibs(1)+1
      call aoib(4)
      do i1=ii,stibs(1)
      stib(i1)=eoa
      enddo
      call trm(4,npmk+1)
      ii=npmk+1
      pmkd(0)=stibs(1)-ii
      pmkl(0)=pmkd(0)-ii
      pmkp(0)=pmkl(0)-ii
      jj=pmkp(0)-llp(0)+1
      call xipht(pmkp(0)+1,stibs(1),-jj)
      call aoib(-jj)
      pmkd(0)=pmkd(0)-jj
      pmkl(0)=pmkl(0)-jj
      pmkp(0)=pmkp(0)-jj
      psize=6+2*npmk
      ii=npmk+1
      apmkd(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      endif
      do i1=1,npmk
      stib(apmkd(0)+i1)=0
      enddo
      if(fpass.eq.1)then
      fpass=0
      goto 50
      endif
      elseif(level.eq.3)then
      if(fpass.eq.1)then
      ii=stibs(1)+1
      call aoib(3)
      do i1=ii,stibs(1)
      stib(i1)=eoa
      enddo
      call trm(3,nvmk+1)
      ii=nvmk+1
      vmkl(0)=stibs(1)-ii
      vmkp(0)=vmkl(0)-ii
      avmkd(0)=vmkp(0)-ii
      call xipht(vmkp(0),stibs(1),-1)
      call aoib(-1)
      vmkl(0)=vmkl(0)-1
      vmkp(0)=vmkp(0)-1
      avmkd(0)=avmkd(0)-1
      else
      nvert=nvert+1
      endif
      do i1=1,nvmk
      stib(avmkd(0)+i1)=0
      enddo
      if(fpass.eq.1)then
      fpass=0
      goto 50
      endif
      endif
      goto 65
  281 continue
      if(idin.eq.3)then
      if(stdw(stcb,id1,id2).eq.0)then
      if(id1.ne.id2)then
      goto 521
      endif
      ii=ichar(stcb(id1:id1))
      if((ii.ne.aplus).and.(ii.ne.aminus))then
      goto 521
      endif
      else
      if(knf.ne.2)then
      auxlin(1:srec)='unknown field in vertex,'
      call messag(1,dbl,nlin,2)
      endif
      knf=0
      call rpi
      level=3
      fpass=1
      lastp=0
      lastkp=0
      goto 50
      endif
      if(knf.ne.0)then
      auxlin(1:srec)='field appeared in previous propagator,'
      call messag(1,dbl,nlin,2)
      endif
      if(fpass.eq.0)then
      ii=ichar(stcb(id1:id1))
      if(ii.eq.aplus)then
      stib(off1+3)=0
      elseif(ii.eq.aminus)then
      stib(off1+3)=1
      else
      goto 521
      endif
      stib(off1+4)=0
      if(off2.ne.0)then
      stib(off2+3)=stib(off1+3)
      stib(off2+4)=0
      endif
      stib(lastp)=off1+1
      stib(off1+1)=eoa
      lastp=off1+1
      if(off2.ne.0)then
      stib(lastp)=off2+1
      stib(off2+1)=eoa
      lastp=off2+1
      endif
      endif
      goto 284
      elseif(idin.eq.4)then
      if(fpass.eq.0)then
      call mstr0(stcb,id1,id2,popt7(0),wopt7(0),ij)
      if(ij.eq.0)then
      goto 521
      endif
      stib(off1+4)=stib(copt7(0)+ij)
      if(off2.ne.0)then
      stib(off2+4)=stib(off1+4)
      endif
      endif
      goto 284
      endif
      if(id1.le.0)then
      goto 521
      elseif(id1.gt.id2)then
      goto 521
      elseif(stdw(stcb,id1,id2).eq.0)then
      auxlin(1:srec)='unacceptable field name,'
      call messag(1,dbl,nlin,2)
      endif
      if(fpass.eq.1)then
      goto 284
      endif
      kk=id2-id1+1
      if(kk.le.0)then
      goto 521
      endif
      newid=1
      if(nphi.gt.0)then
      call mstr1(stcb,id1,id2,llp(0),4,5,ij)
      if(ij.ne.0)then
      newid=0
      endif
      endif
      if(newid.eq.0)then
      if(idin.eq.1)then
      knf=1
      elseif(idin.eq.2)then
      if(knf.eq.1)then
      knf=2
      elseif(ij.ne.nphi)then
      knf=1
      endif
      endif
      endif
      if(knf.ne.0)then
      goto 284
      endif
      if(idin.eq.1)then
      if(nprop.eq.0)then
      call aoib(1)
      stib(stibs(1))=eoa
      llp(0)=stibs(1)
      lastp=stibs(1)
      endif
      off1=stibs(1)
      off2=0
      call aoib(psize)
      do i1=stibs(1)-psize+1,stibs(1)
      stib(i1)=0
      enddo
      stib(off1+5)=id1
      stib(off1+6)=id2-id1+1
      nphi=nphi+1
      elseif(idin.eq.2)then
      newid=1
      i1=stib(off1+5)
      i2=stib(off1+6)
      if(id2-id1+1.eq.i2)then
      if(stcb(id1:id2).eq.stcb(i1:i1-1+i2))then
      newid=0
      endif
      endif
      if(newid.eq.0)then
      stib(off1+2)=0
      goto 284
      endif
      off2=stibs(1)
      call aoib(psize)
      do i1=stibs(1)-psize+1,stibs(1)
      stib(i1)=0
      enddo
      stib(off2+5)=id1
      stib(off2+6)=id2-id1+1
      stib(off1+2)=1
      stib(off2+2)=-1
      nphi=nphi+1
      endif
  284 continue
      id1=0
      goto 287
  282 continue
      if(id1.le.0)then
      goto 521
      elseif(id1.gt.id2)then
      goto 521
      elseif(stdw(stcb,id1,id2).eq.0)then
      auxlin(1:srec)='unacceptable field name,'
      call messag(1,dbl,nlin,2)
      endif
      if(idin.gt.nrho)then
      if(idin.gt.maxdeg)then
      auxlin(1:srec)='vertex degree is too large,'
      call messag(1,dbl,nlin,2)
      endif
      nrho=idin
      endif
      if(fpass.eq.1)then
      goto 286
      endif
      if(idin.eq.1)then
      if(nvert.eq.0)then
      call aoib(1)
      stib(stibs(1))=eoa
      llp(0)=stibs(1)
      lastp=stibs(1)
      endif
      ii=stibs(1)+1
      jj=2+2*nvmk
      call aoib(jj)
      stib(lastp)=ii
      lastp=ii
      stib(lastp)=eoa
      do i1=lastp+1,stibs(1)
      stib(i1)=0
      enddo
      endif
      call mstr0(stcb,id1,id2,namep(0),namel(0),ij)
      if(ij.eq.0)then
      auxlin(1:srec)='unknown field in vertex,'
      call messag(1,dbl,nlin,2)
      endif
      call aoib(1)
      stib(stibs(1))=0
      stib(lastp+1+idin)=ij
      stib(lastp+1)=stib(lastp+1)+1
  286 continue
      id1=0
      goto 287
  283 continue
      if(id1.le.0)then
      goto 521
      elseif(id1.gt.id2)then
      goto 521
      elseif(stdw(stcb,id1,id2).eq.0)then
      goto 521
      endif
      if(level.eq.1)then
      if(kin.gt.1)then
      goto 521
      endif
      if(ngmk.eq.0)then
      call aoib(1)
      stib(stibs(1))=eoa
      llp(0)=stibs(1)
      lastkp=stibs(1)
      else
      call mstr1(stcb,id1,id2,llp(0),1,2,ij)
      if(ij.ne.0)then
      auxlin(1:srec)='multiply defined M-function,'
      call messag(1,dbl,nlin,4-xfile)
      endif
      endif
      ii=stibs(1)+1
      call aoib(4)
      call mstr0(stcb,id1,id2,udkp(0),udkl(0),ij)
      if(ij.eq.0)then
      stib(stibs(1)-2)=id1
      else
      stib(stibs(1)-2)=stib(udkp(0)+ij)
      endif
      stib(stibs(1)-1)=id2-id1+1
      stib(lastkp)=ii
      lastkp=ii
      stib(lastkp)=eoa
      ngmk=ngmk+1
      stib(stibs(1))=0
      elseif(level.eq.2)then
      if(fpass.eq.1)then
      call mstr0(stcb,id1,id2,gmkp(0),gmkl(0),ij)
      if(ij.ne.0)then
      auxlin(1:srec)='multiply defined M-function,'
      call messag(1,dbl,nlin,2)
      endif
      if(npmk.eq.0)then
      call aoib(1)
      stib(stibs(1))=eoa
      llp(0)=stibs(1)
      lastkp=stibs(1)
      else
      call mstr1(stcb,id1,id2,llp(0),1,2,ij)
      if(ij.ne.0)then
      auxlin(1:srec)='multiply defined M-function,'
      call messag(1,dbl,nlin,2)
      endif
      endif
      call mstr0(stcb,id1,id2,udkp(0),udkl(0),ij)
      ii=stibs(1)+1
      call aoib(4)
      if(ij.eq.0)then
      stib(stibs(1)-2)=id1
      stib(stibs(1)-1)=id2-id1+1
      else
      stib(stibs(1)-2)=stib(udkp(0)+ij)
      stib(stibs(1)-1)=stib(udkl(0)+ij)
      endif
      stib(stibs(1))=0
      stib(lastkp)=ii
      lastkp=ii
      stib(lastkp)=eoa
      npmk=npmk+1
      else
      call mstr0(stcb,id1,id2,pmkp(0),pmkl(0),kid)
      if(kid.eq.0)then
      auxlin(1:srec)='unexpected M-function,'
      call messag(1,dbl,nlin,2)
      endif
      if(stib(apmkd(0)+kid).ne.0)then
      auxlin(1:srec)='multiply defined M-function,'
      call messag(1,dbl,nlin,2)
      endif
      endif
      elseif(level.eq.3)then
      if(fpass.eq.1)then
      call mstr0(stcb,id1,id2,gmkp(0),gmkl(0),ij)
      if(ij.ne.0)then
      auxlin(1:srec)='multiply defined M-function,'
      call messag(1,dbl,nlin,2)
      endif
      call mstr0(stcb,id1,id2,pmkp(0),pmkl(0),ij)
      if(ij.ne.0)then
      auxlin(1:srec)='multiply defined M-function,'
      call messag(1,dbl,nlin,2)
      endif
      if(nvmk.eq.0)then
      call aoib(1)
      stib(stibs(1))=eoa
      llp(0)=stibs(1)
      lastkp=stibs(1)
      else
      call mstr1(stcb,id1,id2,llp(0),1,2,ij)
      if(ij.ne.0)then
      auxlin(1:srec)='multiply defined M-function,'
      call messag(1,dbl,nlin,2)
      endif
      endif
      call mstr0(stcb,id1,id2,udkp(0),udkl(0),ij)
      ii=stibs(1)+1
      call aoib(3)
      if(ij.eq.0)then
      stib(stibs(1)-1)=id1
      stib(stibs(1))=id2-id1+1
      else
      stib(stibs(1)-1)=stib(udkp(0)+ij)
      stib(stibs(1))=stib(udkl(0)+ij)
      endif
      stib(lastkp)=ii
      lastkp=ii
      stib(lastkp)=eoa
      nvmk=nvmk+1
      else
      call mstr0(stcb,id1,id2,vmkp(0),vmkl(0),kid)
      if(kid.eq.0)then
      auxlin(1:srec)='unexpected M-function,'
      call messag(1,dbl,nlin,2)
      endif
      if(stib(avmkd(0)+kid).ne.0)then
      auxlin(1:srec)='multiply defined M-function,'
      call messag(1,dbl,nlin,2)
      endif
      endif
      endif
      id1=0
      goto 287
  285 continue
      if(is1.eq.0)then
      goto 521
      elseif(is1.gt.is2)then
      goto 521
      endif
      if(level.eq.1)then
      call aoib(2)
      j1=lastkp+3
      if(vel.eq.0)then
      stib(j1)=1
      else
      stib(j1)=vin+1
      endif
      if(vel.ne.0)then
      auxlin(1:srec)='wrong dimension for M-keyword,'
      call messag(1,dbl,nlin,2)
      endif
      elseif(level.eq.2)then
      if(fpass.eq.1)then
      j1=lastkp+3
      if(vel.eq.0)then
      stib(j1)=1
      else
      stib(j1)=3
      endif
      goto 184
      else
      j1=apmkd(0)+kid
      if(vel.eq.0)then
      stib(j1)=1
      else
      stib(j1)=vin+1
      endif
      if(stib(j1).gt.stib(pmkd(0)+kid))then
      auxlin(1:srec)='wrong dimension for M-function,'
      call messag(1,dbl,nlin,2)
      endif
      endif
      if(vin.gt.2)then
      auxlin(1:srec)='wrong dimension for M-function,'
      call messag(1,dbl,nlin,2)
      endif
      elseif(level.eq.3)then
      if(vel.ne.0)then
      auxlin(1:srec)='wrong dimension for M-function,'
      call messag(1,dbl,nlin,2)
      endif
      if(fpass.eq.1)then
      goto 184
      else
      stib(avmkd(0)+kid)=1
      endif
      endif
      if(ichar(stcb(is1:is1)).eq.squote)then
      kk=stds(stcb,is1,is2,1,1)
      if(kk.lt.0)then
      goto 521
      elseif(kk.gt.0)then
      stcb(is1+1:is1+kk)=stcb(stcbs(1)+1:stcbs(1)+kk)
      endif
      stcb(is1:is1)=char(172)
      stcb(is1+kk+1:is1+kk+1)=char(172)
      is1=is1+1
      ii=-(is2-is1-kk)
      if(ii.lt.0)then
      call cxipht(is2+1,stcbs(1),ii)
      call aocb(ii)
      ix=ix+ii
      top=top+ii
      if(level.eq.2)then
      if(nprop.eq.0)then
      do i1=1,npmk
      if(stib(pmkp(0)+i1).gt.is2)then
      stib(pmkp(0)+i1)=stib(pmkp(0)+i1)+ii
      endif
      enddo
      endif
      elseif(level.eq.3)then
      if(nvert.eq.0)then
      do i1=1,nvmk
      if(stib(vmkp(0)+i1).gt.is2)then
      stib(vmkp(0)+i1)=stib(vmkp(0)+i1)+ii
      endif
      enddo
      endif
      endif
      endif
      else
      if(stdw(stcb,is1,is2).eq.0)then
      if(stdq(stcb,is1,is2).eq.0)then
      goto 521
      endif
      endif
      kk=is2-is1+1
      endif
      if(level.eq.1)then
      stib(stibs(1)-1)=is1
      stib(stibs(1))=kk
      elseif(level.eq.2)then
      jj=5+2*kid
      if(vel.ne.0)then
      if(vin.eq.1)then
      j1=off1+jj
      elseif(vin.eq.2)then
      j1=off2+jj
      endif
      stib(j1)=is1
      stib(j1+1)=kk
      else
      stib(off1+jj)=is1
      stib(off1+jj+1)=kk
      if(off2.ne.0)then
      stib(off2+jj)=is1
      stib(off2+jj+1)=kk
      endif
      endif
      elseif(level.eq.3)then
      jj=lastp+stib(lastp+1)+2*kid
      stib(jj)=is1
      stib(jj+1)=kk
      endif
  184 continue
      is1=0
      goto 287
  521 continue
      auxlin(1:srec)='wrong syntax,'
      call messag(1,dbl,nlin,4-xfile)
  777 continue
      call rvi(llp(0))
      return
      end
      subroutine iki
      implicit integer(a-z)
      save
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( srec=81, ssrec=62 )
      character*(srec) auxlin
      common/z22g/auxlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      common/z38g/gmkp(0:0),gmkl(0:0),gmkd(0:0),gmko(0:0),
     :gmkvp(0:0),gmkvl(0:0)
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      do 100 i=1,nudk
      i1=stib(udkp(0)+i)
      i2=i1-1+stib(udkl(0)+i)
      call mstr0(stcb,i1,i2,gmkp(0),gmkl(0),jj)
      if(jj.ne.0)then
      stib(udkt(0)+i)=1
      stib(udki(0)+i)=jj
      goto 100
      endif
      call mstr0(stcb,i1,i2,pmkp(0),pmkl(0),jj)
      if(jj.ne.0)then
      stib(udkt(0)+i)=2
      stib(udki(0)+i)=jj
      goto 100
      endif
      call mstr0(stcb,i1,i2,vmkp(0),vmkl(0),jj)
      if(jj.eq.0)then
      auxlin(1:srec)='function used in style-file was not found in'
      call messag(1,0,0,2)
      endif
      stib(udkt(0)+i)=3
      stib(udki(0)+i)=jj
  100 continue
      j=udkp(1)
      do 500 i1=1,nudk
      j=stib(j)
      if(stib(udkt(0)+i1).eq.1)then
      do 200 i2=3,7
      if(stib(j+i2).gt.1)then
      auxlin(1:srec)='global M-functions cannot have "dual-" as prefix'
      call messag(1,0,0,4)
      endif
  200 continue
      elseif(stib(udkt(0)+i1).eq.2)then
      if(stib(pmkd(0)+stib(udki(0)+i1)).eq.3)then
      if(stib(j+3).ne.0)then
      auxlin(1:srec)='field M-functions cannot be used outside a loop,'
      call messag(1,0,0,4)
      elseif(stib(j+6).ne.0)then
      auxlin(1:srec)='field M-functions cannot be used in the '//
     :'(strict) vertex loop,'
      call messag(1,0,0,4)
      endif
      else
      if(stib(j+3).ne.0)then
      auxlin(1:srec)='propagator M-functions cannot be used '//
     :'outside a loop,'
      call messag(1,0,0,4)
      elseif(stib(j+6).ne.0)then
      auxlin(1:srec)='propagator M-functions cannot be used in the '//
     :'(strict) vertex loop,'
      call messag(1,0,0,4)
      endif
      do 300 i2=4,7
      if(stib(j+i2).gt.1)then
      auxlin(1:srec)='propagator M-functions cannot be prefixed with '//
     :'"dual-",'
      call messag(1,0,0,4)
      endif
  300 continue
      endif
      elseif(stib(udkt(0)+i1).eq.3)then
      if(stib(j+3).ne.0)then
      auxlin(1:srec)='vertex M-functions cannot be used outside a loop,'
      call messag(1,0,0,4)
      else
      do 400 i2=4,7
      if(i2.ne.5.and.stib(j+i2).gt.1)then
      auxlin(1:srec)='vertex M-functions cannot be prefixed with '//
     :'"dual-" except in the propagator_loop,'
      call messag(1,0,0,4)
      endif
  400 continue
      endif
      endif
  500 continue
      return
      end
      subroutine rsfki
      implicit integer(a-z)
      save
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( eoa=-127, nap=-2047 )
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z36g/udkp(0:2),udkl(0:0),udkt(0:0),udki(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      ii=nudk+1
      udkp(2)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      udkp(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      udkl(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      udkt(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      udki(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      j1=udkp(1)
      do i1=1,nudk
      stib(udkt(0)+i1)=0
      stib(udki(0)+i1)=0
      j1=stib(j1)
      stib(udkp(2)+i1)=j1+1
      stib(udkp(0)+i1)=stib(j1+1)
      stib(udkl(0)+i1)=stib(j1+2)
      enddo
      return
      end
      subroutine rgmki(llpp)
      implicit integer(a-z)
      save
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( eoa=-127, nap=-2047 )
      common/z30g/stib(1:sibuff)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z38g/gmkp(0:0),gmkl(0:0),gmkd(0:0),gmko(0:0),
     :gmkvp(0:0),gmkvl(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      if(ngmk.eq.0)then
      gmkp(0)=nap
      gmkl(0)=nap
      gmkd(0)=nap
      gmko(0)=nap
      gmkvp(0)=nap
      gmkvl(0)=nap
      call aoib(-1)
      goto 90
      endif
      ii=ngmk+1
      gmkp(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      gmkl(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      gmkd(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      gmko(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      j1=0
      tgmkd=0
      jj=llpp
   10 continue
      jj=stib(jj)
      if(jj.ne.eoa)then
      j1=j1+1
      stib(gmkp(0)+j1)=stib(jj+1)
      stib(gmkl(0)+j1)=stib(jj+2)
      stib(gmkd(0)+j1)=stib(jj+3)
      stib(gmko(0)+j1)=tgmkd
      tgmkd=tgmkd+max(1,stib(jj+3)-1)
      goto 10
      endif
      ii=tgmkd+1
      gmkvp(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      gmkvl(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      j1=gmkd(0)
      j2=0
      jj=llpp
   20 continue
      jj=stib(jj)
      if(jj.ne.eoa)then
      j1=j1+1
      j3=jj+3
      do i1=1,max(1,stib(j1)-1)
      j2=j2+1
      j3=j3+1
      stib(gmkvp(0)+j2)=stib(j3)
      j3=j3+1
      stib(gmkvl(0)+j2)=stib(j3)
      if(stib(j3).eq.0)then
      i2=stib(gmkp(0)+j2)
      i3=i2-1+stib(gmkl(0)+j2)
      call wput(4,i2,i3)
      endif
      enddo
      goto 20
      endif
      ii=gmkp(0)-llpp+1
      call xipht(gmkp(0)+1,stibs(1),-ii)
      call aoib(-ii)
      gmkp(0)=gmkp(0)-ii
      gmkl(0)=gmkl(0)-ii
      gmkd(0)=gmkd(0)-ii
      gmko(0)=gmko(0)-ii
      gmkvp(0)=gmkvp(0)-ii
      gmkvl(0)=gmkvl(0)-ii
   90 continue
      return
      end
      subroutine rpi
      implicit integer(a-z)
      save
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( srec=81, ssrec=62 )
      parameter ( eoa=-127, nap=-2047 )
      common/z9g/tpc(0:0)
      common/z11g/nphi,nblok,nprop,npprop
      common/z12g/jflag(1:11),mflag(1:24)
      character*(srec) auxlin
      common/z22g/auxlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z33g/namep(0:0),namel(0:0)
      common/z34g/link(0:0),antiq(0:0),blok(0:0)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z39g/pmkp(0:0),pmkl(0:0),pmkd(0:0),pmkvpp(0:0),pmkvlp(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      common/z42g/pmkr(0:0),pmkvma(0:0),pmkvmi(0:0)
      if(nphi.eq.0)then
      auxlin(1:srec)='no propagators listed in'
      call messag(1,0,0,2)
      endif
      ii=stibs(1)
      psize=6+2*npmk
      call aoib(psize)
      do i1=ii+1,stibs(1)
      stib(i1)=eoa
      enddo
      jj=nphi+1
      call trm(psize,jj)
      blok(0)=stibs(1)-psize*jj
      link(0)=blok(0)+jj
      antiq(0)=link(0)+jj
      tpc(0)=antiq(0)+jj
      namep(0)=tpc(0)+jj
      namel(0)=namep(0)+jj
      if(stib(blok(0)).ne.eoa)then
      stib(blok(0))=eoa
      endif
      npprop=0
      do i1=1,nphi
      j1=stib(tpc(0)+i1)
      j2=stib(link(0)+i1)
      stib(link(0)+i1)=j2+i1
      if(j1.eq.5)then
      mflag(18)=1
      else
      if(j1.eq.1)then
      mflag(19)=1
      endif
      if(j2.ge.0)then
      npprop=npprop+1
      endif
      endif
      enddo
      if(npmk.eq.0)then
      pmkvpp(0)=nap
      pmkvlp(0)=nap
      pmkr(0)=nap
      pmkvmi(0)=nap
      pmkvma(0)=nap
      goto 90
      endif
      ii=npmk+1
      pmkvpp(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      pmkvlp(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      pmkr(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      ii=nphi+1
      jj=namel(0)
      do i1=1,npmk
      jj=jj+ii
      stib(pmkvpp(0)+i1)=jj
      jj=jj+ii
      stib(pmkvlp(0)+i1)=jj
      enddo
      do i1=1,npmk
      ii=stib(pmkvpp(0)+i1)
      jj=stib(pmkvlp(0)+i1)
      ij=0
      do i2=4,2,-2
      if(ij.eq.0)then
      do i3=1,nphi
      j1=stib(ii+i3)
      j2=stib(ii+i3)-1+stib(jj+i3)
      if(j1.gt.j2)then
      ij=1
      goto 40
      elseif(i2.eq.4)then
      if(stdz(stcb,j1,j2).eq.0)then
      goto 40
      endif
      elseif(i2.eq.2)then
      if(stdw(stcb,j1,j2).eq.0)then
      goto 40
      endif
      endif
      enddo
      ij=i2
   40 continue
      endif
      enddo
      if(ij.eq.0)then
      ij=1
      endif
      stib(pmkr(0)+i1)=ij
      enddo
      ii=npmk+1
      pmkvmi(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      pmkvma(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      do i1=1,npmk
      if(stib(pmkr(0)+i1).eq.4)then
      ii=stib(pmkvpp(0)+i1)
      jj=stib(pmkvlp(0)+i1)
      do i2=1,nphi
      j1=stib(ii+i2)
      j2=stib(ii+i2)-1+stib(jj+i2)
      ij=stoz(stcb,j1,j2)
      if(i2.gt.1)then
      if(ij.lt.j3)then
      j3=ij
      elseif(ij.gt.j4)then
      j4=ij
      endif
      else
      j3=ij
      j4=ij
      endif
      enddo
      if(j3.eq.0)then
      if(j4.eq.0)then
      j1=stib(pmkp(0)+i1)
      j2=j1-1+stib(pmkl(0)+i1)
      call wput(6,j1,j2)
      endif
      endif
      else
      j3=0
      j4=0
      endif
      stib(pmkvmi(0)+i1)=j3
      stib(pmkvma(0)+i1)=j4
      enddo
      do i1=1,npmk
      if(stib(pmkr(0)+i1).eq.1)then
      jj=stib(pmkvlp(0)+i1)
      do i2=1,nphi
      if(stib(jj+i2).gt.0)then
      goto 85
      endif
      enddo
      j1=stib(pmkp(0)+i1)
      j2=j1-1+stib(pmkl(0)+i1)
      call wput(5,j1,j2)
      endif
   85 continue
      enddo
   90 continue
      return
      end
      subroutine rvi(llp)
      implicit integer(a-z)
      save
      parameter ( maxleg=11, maxrho=7, maxi=11, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( srec=81, ssrec=62 )
      parameter ( eoa=-127, nap=-2047 )
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1
      common/z3g/nv(0:0),dpntro(0:0),vparto(0:0),vdeg(0:0),nvert
      common/z11g/nphi,nblok,nprop,npprop
      common/z12g/jflag(1:11),mflag(1:24)
      character*(srec) auxlin
      common/z22g/auxlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z34g/link(0:0),antiq(0:0),blok(0:0)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      common/z43g/vmkr(0:0),vmkmao(0:0),vmkmio(0:0)
      integer xtmp(1:maxdeg+1),nrot(1:maxdeg),rotvpo(1:maxdeg)
      if(nvert.eq.0)then
      auxlin(1:srec)='no vertices listed in'
      call messag(1,0,0,2)
      endif
      ii=nrho+1
      nv(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      ii=nvert+1
      vdeg(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      vparto(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      jj=nvmk+1
      vmkvpp(0)=stibs(1)
      call aoib(jj)
      stib(stibs(1))=eoa
      vmkvlp(0)=stibs(1)
      call aoib(jj)
      stib(stibs(1))=eoa
      j1=stibs(1)
      call aoib(2*nvmk*ii)
      do i1=1,nvmk
      stib(vmkvpp(0)+i1)=j1
      j1=j1+ii
      stib(j1)=eoa
      stib(vmkvlp(0)+i1)=j1
      j1=j1+ii
      stib(j1)=eoa
      enddo
      if(j1.ne.stibs(1))then
      auxlin(1:srec)='rvi_4'
      call messag(1,0,0,0)
      endif
      do i1=1,nrho
      nrot(i1)=0
      stib(nv(0)+i1)=0
      enddo
      j1=0
      jj=llp
   10 continue
      jj=stib(jj)
      if(jj.ne.eoa)then
      j1=j1+1
      j2=stib(jj+1)
      stib(vdeg(0)+j1)=j2
      j2=j2+nv(0)
      stib(j2)=stib(j2)+1
      goto 10
      endif
      if((j1.ne.nvert).or.(nrho.lt.3))then
      auxlin(1:srec)='rvi_1'
      call messag(1,0,0,0)
      endif
      j1=0
      jj=llp
   30 continue
      jj=stib(jj)
      if(jj.ne.eoa)then
      j1=j1+1
      j2=jj+stib(vdeg(0)+j1)
      do i1=1,nvmk
      j2=j2+2
      stib(stib(vmkvpp(0)+i1)+j1)=stib(j2)
      stib(stib(vmkvlp(0)+i1)+j1)=stib(j2+1)
      enddo
      goto 30
      endif
      if(j1.ne.nvert)then
      auxlin(1:srec)='rvi_3'
      call messag(1,0,0,0)
      endif
      kk=llp
      if(llp.gt.1)then
      if(stib(llp-1).eq.eoa)then
      kk=llp-1
      endif
      endif
      ij=stib(llp)
      if(stib(llp).ne.eoa)then
      stib(llp)=eoa
      endif
      j1=0
   50 continue
      jj=ij
      if(jj.ne.eoa)then
      j1=j1+1
      ij=stib(ij)
      j2=jj-kk+1
      call xipht(jj+2,jj+1+stib(vdeg(0)+j1),-j2)
      stib(vparto(0)+j1)=kk
      kk=kk+stib(vdeg(0)+j1)
      goto 50
      endif
      if(j1.ne.nvert)then
      auxlin(1:srec)='rvi_5'
      call messag(1,0,0,0)
      endif
      kk=nv(0)-kk-1
      call xipht(nv(0),stibs(1),-kk)
      call aoib(-kk)
      nv(0)=nv(0)-kk
      vdeg(0)=vdeg(0)-kk
      vparto(0)=vparto(0)-kk
      vmkvpp(0)=vmkvpp(0)-kk
      vmkvlp(0)=vmkvlp(0)-kk
      do i1=1,nvmk
      stib(vmkvpp(0)+i1)=stib(vmkvpp(0)+i1)-kk
      stib(vmkvlp(0)+i1)=stib(vmkvlp(0)+i1)-kk
      enddo
      stib(nv(0))=eoa
      jj=nvmk+1
      vmkr(0)=stibs(1)
      call aoib(jj)
      stib(stibs(1))=eoa
      vmks(0)=stibs(1)
      call aoib(jj)
      stib(stibs(1))=eoa
      do i1=vmks(0)+1,vmks(0)+nvmk
      stib(i1)=0
      enddo
      do i1=1,nvmk
      ii=stib(vmkvpp(0)+i1)
      jj=stib(vmkvlp(0)+i1)
      ij=0
      do i2=4,2,-2
      if(ij.eq.0)then
      do i3=1,nvert
      j1=stib(ii+i3)
      j2=stib(ii+i3)-1+stib(jj+i3)
      if(j1.gt.j2)then
      ij=1
      goto 90
      elseif(i2.eq.4)then
      if(stdz(stcb,j1,j2).eq.0)then
      goto 90
      endif
      elseif(i2.eq.2)then
      if(stdw(stcb,j1,j2).eq.0)then
      goto 90
      endif
      endif
      enddo
      ij=i2
   90 continue
      endif
      enddo
      if(ij.eq.0)then
      ij=1
      endif
      stib(vmkr(0)+i1)=ij
      enddo
      do i1=1,nvmk
      if(stib(vmkr(0)+i1).eq.1)then
      jj=stib(vmkvlp(0)+i1)
      do i2=1,nvert
      if(stib(jj+i2).gt.0)then
      goto 23
      endif
      enddo
      j1=stib(vmkp(0)+i1)
      j2=j1-1+stib(vmkl(0)+i1)
      call wput(5,j1,j2)
      elseif(stib(vmkr(0)+i1).eq.4)then
      ii=stib(vmkvpp(0)+i1)
      jj=stib(vmkvlp(0)+i1)
      do i2=1,nvert
      j1=stib(ii+i2)
      j2=j1-1+stib(jj+i2)
      if(nullz(stcb,j1,j2).eq.0)then
      goto 23
      endif
      enddo
      j1=stib(vmkp(0)+i1)
      j2=j1-1+stib(vmkl(0)+i1)
      call wput(6,j1,j2)
      endif
   23 continue
      enddo
      nvrot=0
      do 403 i1=1,nrho
      if(stib(nv(0)+i1).gt.0)then
      rotvpo(i1)=stibs(1)
      do i2=1,nvert
      if(i1.eq.stib(vdeg(0)+i2))then
      j1=stib(vparto(0)+i2)
      do i3=1,i1
      xtmp(i3)=stib(j1+i3)
      enddo
      do i3=1,i1-1
      do i4=i3+1,i1
      if(xtmp(i3).gt.xtmp(i4))then
      ii=xtmp(i3)
      xtmp(i3)=xtmp(i4)
      xtmp(i4)=ii
      endif
      enddo
      enddo
      goto 830
  805 continue
      do i3=i1-1,1,-1
      if(xtmp(i3).lt.xtmp(i3+1))then
      j1=i3
      goto 810
      endif
      enddo
      goto 401
  810 continue
      do i3=j1+2,i1
      if(xtmp(i3).le.xtmp(j1))then
      j2=i3-1
      goto 820
      endif
      enddo
      j2=i1
  820 continue
      ii=xtmp(j1)
      xtmp(j1)=xtmp(j2)
      xtmp(j2)=ii
      jj=(i1-j1)/2
      j1=j1+1
      do i3=0,jj-1
      ii=xtmp(j1+i3)
      xtmp(j1+i3)=xtmp(i1-i3)
      xtmp(i1-i3)=ii
      enddo
  830 continue
      j1=stibs(1)+1
      call aoib(i1+1)
      stib(j1)=i2
      do i3=1,i1
      stib(j1+i3)=xtmp(i3)
      enddo
      nvrot=nvrot+1
      nrot(i1)=nrot(i1)+1
      goto 805
      endif
  401 continue
      enddo
      j1=stibs(1)+1
      call aoib(i1+1)
      stib(j1)=0
      do i2=1,i1
      stib(j1+i2)=nphi+1
      enddo
      endif
  403 continue
      call aoib(1)
      stib(stibs(1))=eoa
      do 790 i1=1,nrho
      if((stib(nv(0)+i1).gt.0).and.(nrot(i1).gt.1))then
      h2=nrot(i1)
      h1=(h2/2)+1
  720 continue
      if(h1.gt.1)then
      h1=h1-1
      ii=rotvpo(i1)+(i1+1)*(h1-1)
      do i2=1,i1+1
      xtmp(i2)=stib(ii+i2)
      enddo
      else
      ii=rotvpo(i1)+(i1+1)*(h2-1)
      do i2=1,i1+1
      xtmp(i2)=stib(ii+i2)
      enddo
      do i2=1,i1+1
      stib(ii+i2)=stib(rotvpo(i1)+i2)
      enddo
      h2=h2-1
      if(h2.eq.1)then
      do i2=1,i1+1
      stib(rotvpo(i1)+i2)=xtmp(i2)
      enddo
      goto 790
      endif
      endif
      hj=h1
  740 continue
      hi=hj
      hj=hi+hi
      if(hj.le.h2)then
      if(hj.lt.h2)then
      kk=0
      jj=rotvpo(i1)+(i1+1)*(hj-1)
      jk=jj+i1+1
      j1=1
  745 continue
      if(j1.le.i1)then
      j1=j1+1
      if(stib(jj+j1).lt.stib(jk+j1))then
      kk=-1
      elseif(stib(jj+j1).gt.stib(jk+j1))then
      kk=1
      else
      goto 745
      endif
      endif
      if(kk.eq.-1)then
      hj=hj+1
      elseif(kk.eq.0)then
      mflag(17)=1
      endif
      endif
      kk=0
      jj=rotvpo(i1)+(i1+1)*(hj-1)
      j1=1
  755 continue
      if(j1.le.i1)then
      j1=j1+1
      if(xtmp(j1).lt.stib(jj+j1))then
      kk=-1
      elseif(xtmp(j1).gt.stib(jj+j1))then
      kk=1
      else
      goto 755
      endif
      endif
      if(kk.eq.-1)then
      ii=rotvpo(i1)+(i1+1)*(hi-1)
      do i2=1,i1+1
      stib(ii+i2)=stib(jj+i2)
      enddo
      goto 740
      elseif(kk.eq.0)then
      mflag(17)=1
      endif
      endif
      ii=rotvpo(i1)+(i1+1)*(hi-1)
      do i2=1,i1+1
      stib(ii+i2)=xtmp(i2)
      enddo
      goto 720
      endif
  790 continue
      if(mflag(17).ne.0)then
      call wput(3,0,0)
      endif
      do i1=1,nrho
      if(nrot(i1).gt.0)then
      ii=rotvpo(i1)
      do i2=1,nrot(i1)
      j1=ii+2
      j2=j1+i1+1
      jj=i1
   84 continue
      if(jj.gt.0)then
      kk=stib(j2)-stib(j1)
      if(kk.eq.0)then
      jj=jj-1
      j1=j1+1
      j2=j2+1
      goto 84
      elseif(kk.lt.0)then
      auxlin(1:srec)='rvi_2'
      call messag(1,0,0,0)
      endif
      endif
      ii=ii+i1+1
      enddo
      endif
      enddo
      ii=nrho+1
      dpntro(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      do i1=1,nrho
      if(nrot(i1).gt.0)then
      stib(dpntro(0)+i1)=stibs(1)
      call aoib(nphi+1)
      stib(stibs(1))=nap
      ii=rotvpo(i1)+2
      jj=0
      do i2=1,nrot(i1)+1
      if(stib(ii).ne.jj)then
      if(stib(ii).le.nphi)then
      kk=stib(ii)
      else
      kk=nphi
      endif
      do i3=jj+1,kk
      stib(stib(dpntro(0)+i1)+i3)=ii-1
      enddo
      jj=stib(ii)
      endif
      ii=ii+i1+1
      enddo
      else
      stib(dpntro(0)+i1)=nap
      endif
      enddo
      stib(stibs(1))=eoa
      bloka=stibs(1)
      ii=nphi
      call vaoib(ii)
      do i1=blok(0)+1,blok(0)+nphi
      stib(i1)=0
      enddo
      nblok=0
      j1=0
      j2=0
   31 continue
      if(j1.lt.nphi)then
      nblok=nblok+1
      j3=nblok
   32 continue
      if(stib(blok(0)+j3).gt.0)then
      j3=j3+1
      goto 32
      endif
      j1=j1+1
      j2=j1
      stib(bloka+j1)=j3
      stib(blok(0)+j3)=nblok
      ii=stib(link(0)+j3)
      jj=blok(0)+ii
      if(stib(jj).eq.0)then
      j1=j1+1
      stib(bloka+j1)=ii
      stib(jj)=nblok
      endif
   33 continue
      if(j2.le.j1)then
      jj=stib(bloka+j2)
      do i1=2,nrho
      if(stib(dpntro(0)+i1).gt.0)then
      j3=stib(stib(dpntro(0)+i1)+jj)
   36 continue
      if(stib(j3+1).eq.jj)then
      ii=stib(j3+2)
      if(stib(blok(0)+ii).eq.0)then
      j1=j1+1
      stib(bloka+j1)=ii
      stib(blok(0)+ii)=nblok
      endif
      ii=stib(link(0)+ii)
      if(stib(blok(0)+ii).eq.0)then
      j1=j1+1
      stib(bloka+j1)=ii
      stib(blok(0)+ii)=nblok
      endif
      j3=j3+i1+1
      goto 36
      endif
      endif
      enddo
      j2=j2+1
      goto 33
      endif
      goto 31
      endif
      if(nblok.gt.1)then
      call wput(1,0,0)
      endif
      j1=0
   34 continue
      if(j1.lt.nphi)then
      j1=j1+1
      do i1=1,nrho
      if(stib(nv(0)+i1).gt.0)then
      if(stib(stib(stib(dpntro(0)+i1)+j1)+1).eq.j1)then
      goto 34
      endif
      endif
      enddo
      call wput(2,0,0)
      endif
      return
      end
      subroutine vsig
      implicit integer(a-z)
      save
      parameter ( maxleg=11, maxrho=7, maxi=11, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( srec=81, ssrec=62 )
      parameter ( eoa=-127, nap=-2047 )
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1
      common/z3g/nv(0:0),dpntro(0:0),vparto(0:0),vdeg(0:0),nvert
      common/z11g/nphi,nblok,nprop,npprop
      common/z12g/jflag(1:11),mflag(1:24)
      character*(srec) auxlin
      common/z22g/auxlin
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z40g/vmkp(0:0),vmkl(0:0),vmkvpp(0:0),vmkvlp(0:0),vmks(0:0)
      common/z41g/nudk,ngmk,npmk,nvmk,tgmkd,tpmkd
      common/z43g/vmkr(0:0),vmkmao(0:0),vmkmio(0:0)
      if(mflag(20).eq.0)then
      vmkmio(0)=nap
      vmkmao(0)=nap
      goto 91
      endif
      ii=nvmk+1
      vmkmio(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      vmkmao(0)=stibs(1)
      call aoib(ii)
      stib(stibs(1))=eoa
      do i1=1,nvmk
      if(stib(vmks(0)+i1).ne.1)then
      stib(vmkmio(0)+i1)=nap
      stib(vmkmao(0)+i1)=nap
      else
      ii=nrho
      p1=stibs(1)
      stib(vmkmio(0)+i1)=p1
      call aoib(ii)
      p2=stibs(1)
      stib(vmkmao(0)+i1)=p2
      call aoib(ii)
      call vaoib(ii)
      do i2=stibs(1)+1,stibs(1)+ii
      stib(i2)=0
      enddo
      ii=stib(vmkvpp(0)+i1)
      jj=stib(vmkvlp(0)+i1)
      do i2=1,nvert
      j1=stib(ii+i2)
      j2=j1-1+stib(jj+i2)
      ij=stoz(stcb,j1,j2)
      j1=stib(vdeg(0)+i2)
      j2=stibs(1)+j1
      if(stib(j2).ne.0)then
      if(ij.gt.stib(p2+j1))then
      stib(p2+j1)=ij
      elseif(ij.lt.stib(p1+j1))then
      stib(p1+j1)=ij
      endif
      else
      stib(j2)=1
      stib(p1+j1)=ij
      stib(p2+j1)=ij
      endif
      enddo
      endif
      enddo
      if(stib(stibs(1)).ne.eoa)then
      call aoib(1)
      stib(stibs(1))=eoa
      endif
   91 continue
      return
      end
      subroutine gen10(cntr)
      implicit integer(a-z)
      save
      parameter ( maxleg=11, maxrho=7, maxi=11, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( ipar1= maxleg*maxleg-maxleg )
      parameter ( ipar2= ipar1/2 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( srec=81, ssrec=62 )
      common/z1in/leg(1:maxleg),nleg,nloop,incom
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1
      common/z4g/n,nli
      common/z6g/p1(1:maxleg),invp1(1:maxleg)
      common/z7g/lmap(1:maxn,1:maxdeg),vmap(1:maxn,1:maxdeg),
     :pmap(1:maxn,1:maxdeg),vlis(1:maxn),invlis(1:maxn)
      common/z8g/degree(1:maxn),xn(1:maxn)
      common/z10g/p1l(1:ipar2),p1r(1:ipar2),ns1
      common/z12g/jflag(1:11),mflag(1:24)
      common/z16g/rdeg(1:maxn),amap(1:maxn,1:maxdeg)
      common/z18g/eg(1:maxn,1:maxn),flow(1:maxli,0:maxleg+maxrho)
      common/z20g/tftyp(0:0),tfnarg(0:0),tfa(0:0),tfb(0:0),tfc(0:0),
     :tfo(0:0),tf2(0:0),ntf
      character*(srec) auxlin
      common/z22g/auxlin
      common/z30g/stib(1:sibuff)
      common/z50g/tfta(0:0),tftb(0:0),tftic(0:0),ntft
      integer vaux(1:maxli),xli(1:maxleg),xlj(1:maxn)
      if(cntr.eq.0)then
      if(mflag(11).eq.0)then
      goto 29
      else
      goto 27
      endif
      endif
      cntr21=1
   27 continue
      call gen21(cntr21)
      if(cntr21.eq.0)then
      cntr=0
      goto 90
      endif
      do i1=1,n
      vaux(i1)=xn(i1)
      enddo
      do i1=1,rho(1)
      vlis(i1)=i1
      invlis(i1)=i1
      enddo
      do i1=rhop1,n
      invlis(i1)=0
      enddo
      aux=rho(1)
      jj=rhop1
   10 continue
      if(aux.lt.n)then
   20 continue
      if(invlis(jj).ne.0)then
      jj=jj+1
      goto 20
      endif
      ii=jj
      do i1=ii+1,n
      if(invlis(i1).eq.0)then
      if(vaux(i1).gt.vaux(ii))then
      ii=i1
      endif
      endif
      enddo
      if(vaux(ii).eq.0)then
      auxlin(1:srec)='gen10_1'
      call messag(1,0,0,0)
      endif
      aux=aux+1
      vlis(aux)=ii
      invlis(ii)=aux
      rdeg(ii)=vaux(ii)
      do i1=rhop1,n
      vaux(i1)=vaux(i1)+g(i1,ii)
      enddo
      goto 10
      endif
      do i1=1,n
      vaux(i1)=0
      enddo
      do i1=1,n
      ii=vlis(i1)
      kk=0
      j1=1
   30 continue
      if(kk.lt.degree(ii))then
      jj=vlis(j1)
      aux=1
      do i3=1,g(ii,jj)
      kk=kk+1
      vmap(ii,kk)=jj
      vaux(jj)=vaux(jj)+1
      if(ii.ne.jj)then
      lmap(ii,kk)=vaux(jj)
      else
      lmap(ii,kk)=vaux(jj)+aux
      aux=-aux
      endif
      enddo
      j1=j1+1
      goto 30
      endif
      if(kk.ne.degree(ii))then
      auxlin(1:srec)='gen10_2'
      call messag(1,0,0,0)
      endif
      enddo
      do i1=rhop1,n
      do i2=1,degree(i1)-1
      if(invlis(vmap(i1,i2)).gt.invlis(vmap(i1,i2+1)))then
      auxlin(1:srec)='gen10_3'
      call messag(1,0,0,0)
      endif
      enddo
      enddo
      do i1=1,n
      do i2=1,degree(i1)
      j1=vmap(i1,i2)
      j2=lmap(i1,i2)
      if(vmap(j1,j2).ne.i1)then
      auxlin(1:srec)='gen10_4'
      call messag(1,0,0,0)
      endif
      if(lmap(j1,j2).ne.i2)then
      auxlin(1:srec)='gen10_5'
      call messag(1,0,0,0)
      endif
      enddo
      enddo
      if(jflag(1)+jflag(10).ne.0)then
      do i1=rhop1,n
      do i2=1,degree(i1)
      j3=vmap(i1,i2)
      if(j3.gt.rho(1))then
      j1=min(i1,j3)
      j2=max(i1,j3)
      ii=eg(j1,j2)
      jj=g(i1,j3)-1
      if(jj.eq.0)then
      amap(i1,i2)=ii
      elseif(i1.eq.j3)then
      amap(i1,i2)=ii+(i2-1-rdeg(i1))/2
      else
      kk=0
      do i3=i2-1,1,-1
      if(vmap(i1,i3).eq.j3)then
      kk=kk+1
      else
      goto 37
      endif
      enddo
   37 continue
      amap(i1,i2)=ii+kk
      endif
      endif
      enddo
      enddo
      endif
      do i1=1,rho(1)
      p1(i1)=i1
      enddo
      goto 121
   29 continue
      do i1=rho(1)-1,1,-1
      if(p1(i1).lt.p1(i1+1))then
      j1=i1
      goto 38
      endif
      enddo
      goto 27
   38 continue
      do i1=j1+2,rho(1)
      if(p1(i1).lt.p1(j1))then
      j2=i1-1
      goto 61
      endif
      enddo
      j2=rho(1)
   61 continue
      ii=p1(j1)
      p1(j1)=p1(j2)
      p1(j2)=ii
      j1=j1+1
      j2=rho(1)
   33 continue
      if(j1.lt.j2)then
      ii=p1(j1)
      p1(j1)=p1(j2)
      p1(j2)=ii
      j1=j1+1
      j2=j2-1
      goto 33
      endif
  121 continue
      do i1=1,ns1
      if(p1(p1r(i1)).lt.p1(p1l(i1)))then
      do i2=p1r(i1)+1,rho(1)-1
      do i3=i2+1,rho(1)
      if(p1(i2).lt.p1(i3))then
      ii=p1(i2)
      p1(i2)=p1(i3)
      p1(i3)=ii
      endif
      enddo
      enddo
      goto 29
      endif
      enddo
      do i1=1,rho(1)
      invp1(p1(i1))=i1
      enddo
      if(mflag(22).ne.0)then
      ii=stib(tftic(0)+11)
      do i1=stib(tfta(0)+ii),stib(tftb(0)+ii)
      i3=stib(tf2(0)+i1)
      styp=stib(tftyp(0)+i3)
      j1=stib(tfnarg(0)+i3)
      j2=stib(tfo(0)+i3)
      do i2=rhop1,n
      xlj(i2)=0
      enddo
      ii=0
      do i2=1,j1
      j3=vmap(invp1(stib(j2+i2)),1)
      if(xlj(j3).eq.0)then
      xlj(j3)=1
      ii=ii+1
      endif
      enddo
      if((ii.lt.stib(tfa(0)+i3)).or.(ii.gt.stib(tfb(0)+i3)))then
      if(styp.gt.0)then
      goto 29
      else
      goto 55
      endif
      endif
      if(stib(tfc(0)+i3).eq.0)then
      if(styp.gt.0)then
      goto 55
      else
      goto 29
      endif
      endif
      do i2=1,rho(1)
      xli(i2)=0
      enddo
      do i2=1,j1
      xli(invp1(stib(j2+i2)))=1
      enddo
      do i2=1,rho(1)
      if(xli(i2).eq.0)then
      if(xlj(vmap(i2,1)).ne.0)then
      if(styp.gt.0)then
      goto 29
      else
      goto 55
      endif
      endif
      endif
      enddo
      if(styp.lt.0)then
      goto 29
      endif
   55 continue
      enddo
      endif
      if(mflag(23).ne.0)then
      ii=stib(tftic(0)+12)
      do i1=stib(tfta(0)+ii),stib(tftb(0)+ii)
      i3=stib(tf2(0)+i1)
      styp=stib(tftyp(0)+i3)
      j1=stib(tfnarg(0)+i3)
      j2=stib(tfo(0)+i3)
      do i2=1,rho(1)
      xli(i2)=0
      enddo
      do i2=1,j1
      xli(invp1(stib(j2+i2)))=1
      enddo
      jj=0
      do i2=rhop1,nli
      if(flow(i2,0).eq.2)then
      kk=0
      do i3=1,rho(1)
      if(flow(i2,i3).eq.0)then
      if(xli(i3).ne.0)then
      kk=kk+1
      endif
      else
      if(xli(i3).eq.0)then
      kk=kk+1
      endif
      endif
      enddo
      if(kk.gt.0)then
      kk=kk-rho(1)
      endif
      if(kk.eq.0)then
      if(styp.gt.0)then
      goto 57
      else
      goto 29
      endif
      endif
      endif
      enddo
      if(styp.gt.0)then
      goto 29
      endif
   57 continue
      enddo
      endif
      cntr=1
      if(jflag(1).eq.0)then
      goto 90
      endif
      do i1=1,rho(1)
      amap(i1,1)=i1
      amap(vmap(i1,1),lmap(i1,1))=i1
      enddo
      do i1=1,nli
      vaux(i1)=-2
      enddo
      do i1=1,n
      do i2=1,degree(i1)
      ii=amap(i1,i2)
      if(ii.gt.0)then
      if(ii.le.nli)then
      vaux(ii)=vaux(ii)+1
      endif
      endif
      enddo
      enddo
      do i1=1,nli
      if(vaux(i1).ne.0)then
      auxlin(1:srec)='gen10_6'
      call messag(1,0,0,0)
      endif
      enddo
   90 continue
      return
      end
      subroutine gen21(cntr21)
      implicit integer(a-z)
      save
      parameter ( maxleg=11, maxrho=7, maxi=11, maxdeg=6 )
      parameter ( maxn=maxi+maxi-2, maxli=maxn+maxrho )
      parameter ( ipar1= maxleg*maxleg-maxleg )
      parameter ( ipar2= ipar1/2 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( srec=81, ssrec=62 )
      parameter ( eoa=-127, nap=-2047 )
      common/z1g/g(1:maxn,1:maxn),rho(1:maxdeg),nrho,rhop1
      common/z4g/n,nli
      common/z5g/psym(0:0),psyms,nsym
      common/z7g/lmap(1:maxn,1:maxdeg),vmap(1:maxn,1:maxdeg),
     :pmap(1:maxn,1:maxdeg),vlis(1:maxn),invlis(1:maxn)
      common/z8g/degree(1:maxn),xn(1:maxn)
      common/z10g/p1l(1:ipar2),p1r(1:ipar2),ns1
      common/z12g/jflag(1:11),mflag(1:24)
      common/z14g/zcho(0:maxli),zbri(0:maxli),zpro(0:maxli),
     :rbri(0:maxli),sbri(0:maxli)
      common/z17g/xtail(1:maxn),xhead(1:maxn),ntadp
      common/z18g/eg(1:maxn,1:maxn),flow(1:maxli,0:maxleg+maxrho)
      character*(srec) auxlin
      common/z22g/auxlin
      common/z30g/stib(1:sibuff)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      integer xc(1:maxdeg)
      integer xl(1:maxn),xt(1:maxn)
      integer bound(1:maxdeg),degr(1:maxdeg),xs(1:maxn)
      integer xg(1:maxn,1:maxn),ds(1:maxn,1:maxn)
      integer p1s(1:maxleg,1:maxleg)
      integer aa(1:maxn),bb(1:maxn)
      integer uset(0:maxn),xset(1:maxn),xp(1:maxn),a1(1:maxn)
      integer str(1:maxn),dta(1:maxn),lps(1:maxn)
      integer head(1:maxli),tail(1:maxli),intree(1:maxli)
      integer emul(1:maxdeg,1:maxli),nemul(1:maxdeg)
      if(abs(cntr21).ne.1)then
      goto 97
      elseif(nrho.le.0)then
      goto 97
      endif
      do i1=1,nrho
      if((rho(i1).lt.0).or.(rho(i1).gt.maxn))then
      goto 97
      endif
      enddo
      if(cntr21.eq.-1)then
      goto 05
      endif
      n=0
      nli=0
      do i1=1,nrho
      nli=nli+i1*rho(i1)
      do i2=1,rho(i1)
      n=n+1
      degree(n)=i1
      enddo
      enddo
      if((n.le.0).or.(n.gt.maxn))then
      goto 97
      endif
      ii=nli/2
      if(nli.ne.ii+ii)then
      auxlin(1:srec)='gen21_1'
      call messag(1,0,0,0)
      endif
      nli=ii
      loop=nli-n+1
      if(loop.lt.0)then
      goto 97
      endif
      j1=0
      if(mflag(9).gt.0)then
      j1=1
      endif
      do i1=rhop1,n
      if(j1.ne.0)then
      str(i1)=1
      elseif(degree(i1).ne.degree(n))then
      str(i1)=min(degree(i1),loop+1)
      elseif(n.gt.2)then
      str(i1)=min(degree(i1)-1,loop+1)
      else
      str(i1)=min(degree(i1),loop+1)
      endif
      enddo
      j1=1
      xl(1)=1
      do i1=2,n
      if(degree(i1).ne.degree(i1-1))then
      xt(j1)=i1-1
      j1=j1+1
      xl(j1)=i1
      endif
      enddo
      xt(j1)=n
      cntr21=-1
      do i1=1,n
      xn(i1)=0
      enddo
      if(rho(1).eq.0)then
      xc(1)=0
      goto 59
      endif
      ubound=0
      n1=1
      do i1=2,nrho
      if(rho(i1).gt.0)then
      n1=n1+1
      degr(n1)=i1
      bound(n1)=i1*rho(i1)
      ubound=ubound+bound(n1)
      endif
      enddo
      xc(1)=max(0,rho(1)-ubound)
   61 continue
      ii=rho(1)-xc(1)
      do i1=n1,2,-1
      xc(i1)=min(ii,bound(i1))
      ii=ii-xc(i1)
      enddo
   43 continue
      do i1=n1,2,-1
      ii=xc(i1)
      do i2=xl(i1),xt(i1)
      xn(i2)=min(ii,degr(i1))
      ii=ii-xn(i2)
      xs(i2)=ii
      enddo
      enddo
      goto 12
   76 continue
      do i1=n1,2,-1
      do i2=xt(i1)-1,xl(i1),-1
      if(xn(i2).gt.1)then
      j1=i2
      goto 89
      endif
      enddo
      goto 82
   89 continue
      xn(j1)=xn(j1)-1
      xs(j1)=xs(j1)+1
      do i2=j1+1,xt(i1)
      xn(i2)=min(xn(j1),xs(i2-1))
      xs(i2)=xs(i2-1)-xn(i2)
      enddo
      if(xs(xt(i1)).gt.0)then
      do i2=j1-1,xl(i1),-1
      if(xn(i2).gt.1)then
      j1=i2
      goto 89
      endif
      enddo
      goto 82
      endif
      do i2=i1+1,n1
      ii=xc(i2)
      do i3=xl(i2),xt(i2)
      xn(i3)=min(ii,degr(i2))
      ii=ii-xn(i3)
      xs(i3)=ii
      enddo
      enddo
      goto 12
   82 continue
      enddo
      do i1=n1,3,-1
      if(xc(i1).gt.0)then
      j1=i1
      goto 45
      endif
      enddo
      goto 85
   45 continue
      do i1=j1-1,2,-1
      if(xc(i1).lt.bound(i1))then
      j1=i1
      goto 55
      endif
      enddo
      goto 85
   55 continue
      xc(j1)=xc(j1)+1
      ii=-1
      do i1=j1+1,n1
      ii=ii+xc(i1)
      enddo
      do i1=n1,j1+1,-1
      xc(i1)=min(ii,bound(i1))
      ii=ii-xc(i1)
      enddo
      if(ii.ne.0)then
      auxlin(1:srec)='gen21_2'
      call messag(1,0,0,0)
      endif
      goto 43
   85 continue
      xc(1)=xc(1)+2
      if(xc(1).le.rho(1))then
      goto 61
      endif
      goto 80
   12 continue
      if(rho(1).lt.n)then
      if(xc(1).gt.0)then
      goto 80
      endif
      endif
      if(n.gt.rhop1)then
      j1=0
      if(mflag(1).gt.0)then
      j1=1
      endif
      do i1=rhop1,n
      if(xn(i1)+j1.ge.degree(i1))then
      goto 76
      endif
      enddo
      endif
      ii=1
      if(loop.eq.0)then
      ii=0
      elseif((mflag(8).gt.0).and.(rhop1.ne.2))then
      ii=0
      elseif(mflag(9).gt.0)then
      ii=0
      elseif(mflag(10).gt.0)then
      ii=0
      else
      endif
      if(ii.eq.0)then
      do i1=rhop1,n
      dta(i1)=0
      enddo
      else
      do i1=rhop1,n
      ii=0
      if(n.gt.rhop1)then
      if(mflag(1).gt.0)then
      ii=2
      elseif((mflag(2).gt.0).and.(xn(i1).eq.0))then
      ii=2
      else
      ii=1
      endif
      endif
      j=max(0,min(degree(i1)-xn(i1)-ii,loop+loop))/2
      dta(i1)=j+j
      enddo
      endif
   59 continue
      do i1=1,n
      do i2=i1+1,n
      xg(i1,i2)=0
      enddo
      enddo
      do i1=2,xc(1),2
      xg(i1-1,i1)=1
      enddo
      limin=rhop1
      limax=max(limin,n-1)
      lin=limin
      j1=xc(1)
      do i1=limin,n
      if(xn(i1).gt.0)then
      a1(i1)=j1+1
      do i2=1,xn(i1)
      j1=j1+1
      xg(j1,i1)=1
      vmap(j1,1)=i1
      enddo
      else
      a1(i1)=0
      endif
      enddo
      dsum=-1
   70 continue
      dsum=dsum+1
      if(dsum.gt.loop)then
      goto 76
      endif
      ii=2*dsum
      do i1=n,rhop1,-1
      jj=min(ii,dta(i1))
      xg(i1,i1)=jj
      ii=ii-jj
      enddo
      if(ii.eq.0)then
      goto 19
      else
      goto 76
      endif
   32 continue
      do i1=n,rhop1,-1
      if(xg(i1,i1).gt.0)then
      j1=i1
      goto 33
      endif
      enddo
      goto 70
   33 continue
      do i1=j1-1,rhop1,-1
      if(xg(i1,i1).lt.dta(i1))then
      xg(i1,i1)=xg(i1,i1)+2
      j2=i1
      goto 64
      endif
      enddo
      goto 70
   64 continue
      ii=-2
      do i1=j2+1,j1
      ii=ii+xg(i1,i1)
      enddo
      do i1=n,j2+1,-1
      jj=min(ii,dta(i1))
      xg(i1,i1)=jj
      ii=ii-jj
      enddo
      if(ii.ne.0)then
      auxlin(1:srec)='gen21_3'
      call messag(1,0,0,0)
      endif
   19 continue
      do i1=limin,n
      ds(limin,i1)=degree(i1)-xn(i1)-xg(i1,i1)
      enddo
      uset(0)=0
      xset(1)=1
      jj=1
      do i1=2,n
      if(degree(i1-1).ne.degree(i1))then
      uset(jj)=i1-1
      jj=jj+1
      elseif(xn(i1-1).ne.xn(i1))then
      uset(jj)=i1-1
      jj=jj+1
      else
      ii=xg(i1-1,i1-1)-xg(i1,i1)
      if(ii.gt.0)then
      uset(jj)=i1-1
      jj=jj+1
      elseif(ii.lt.0)then
      goto 32
      endif
      endif
      xset(i1)=jj
      enddo
      uset(jj)=n
      lps(lin)=loop-dsum+1
   10 continue
      ii=ds(lin,lin)
      bond=min(str(lin),lps(lin))
      do i1=n,lin+1,-1
      jj=min(ii,bond,ds(lin,i1))
      xg(lin,i1)=jj
      ii=ii-jj
      enddo
      if(ii.gt.0)then
      goto 15
      endif
      goto 28
   05 continue
      lin=limax
      goto 17
   15 continue
      if(lin.eq.limin)then
      goto 32
      endif
      lin=lin-1
   17 continue
      if((lin.lt.limin).or.(lin.gt.n))then
      auxlin(1:srec)='gen21_5'
      call messag(1,0,0,0)
      endif
      do col=n,lin+1,-1
      aux=xg(lin,col)-1
      if(aux.ge.0)then
      goto 23
      endif
      enddo
      goto 15
   23 continue
      bond=min(str(lin),lps(lin))
      do i1=col-1,lin+1,-1
      if(min(ds(lin,i1),bond).gt.xg(lin,i1))then
      xg(lin,i1)=xg(lin,i1)+1
      do i2=n,i1+1,-1
      xg(lin,i2)=min(aux,bond,ds(lin,i2))
      aux=aux-xg(lin,i2)
      enddo
      goto 28
      endif
      aux=aux+xg(lin,i1)
      enddo
      goto 15
   38 continue
      aux=-1
      do i1=col,n
      aux=aux+xg(lin,i1)
      enddo
      goto 23
   28 continue
      if(lin.eq.n)then
      goto 200
      endif
      msum=0
      do i1=lin+1,n
      ii=xg(lin,i1)-1
      if(ii.gt.0)then
      msum=msum+ii
      endif
      enddo
      if(msum.ge.lps(lin))then
      goto 17
      endif
      if(lin.gt.limin)then
      if(xset(lin).eq.xset(lin-1))then
      do i1=limin,lin-2
      ii=xg(i1,lin-1)-xg(i1,lin)
      if(ii.gt.0)then
      goto 35
      elseif(ii.lt.0)then
      auxlin(1:srec)='gen21_7'
      call messag(1,0,0,0)
      endif
      enddo
      do col=lin+1,n
      ii=xg(lin-1,col)-xg(lin,col)
      if(ii.lt.0)then
      goto 38
      elseif(ii.gt.0)then
      goto 35
      endif
      enddo
      endif
      endif
   35 continue
      do col=lin+2,n
      if(xset(col).eq.xset(col-1))then
      do i1=limin,lin
      ii=xg(i1,col-1)-xg(i1,col)
      if(ii.lt.0)then
      goto 38
      elseif(ii.gt.0)then
      goto 24
      endif
      enddo
      endif
   24 continue
      enddo
      do i1=lin+1,n
      ds(lin+1,i1)=ds(lin,i1)-xg(lin,i1)
      if(ds(lin+1,i1).lt.0)then
      auxlin(1:srec)='gen21_8'
      call messag(1,0,0,0)
      endif
      enddo
      lin=lin+1
      lps(lin)=lps(lin-1)-msum
      goto 10
  200 continue
      j1=rhop1
      j2=n-1
      lin=n
      do i1=j1,j2
      aa(i1)=0
      enddo
      aa(n)=1
      i1=n-rho(1)
      bb(i1)=lin
      kk=i1-1
   21 continue
      if(kk.gt.0)then
      ii=bb(i1)
      i1=i1-1
      do i2=j2,ii+1,-1
      if(aa(i2).eq.0)then
      if(xg(ii,i2).ne.0)then
      bb(kk)=i2
      kk=kk-1
      aa(i2)=-aa(ii)
      endif
      endif
      enddo
      if(aa(j2).ne.0)then
      j2=j2-1
      endif
      do i2=j1,ii-1
      if(aa(i2).eq.0)then
      if(xg(i2,ii).ne.0)then
      bb(kk)=i2
      kk=kk-1
      aa(i2)=-aa(ii)
      endif
      endif
      enddo
      if(aa(j1).ne.0)then
      j1=j1+1
      endif
      if(i1.eq.kk)then
      do i2=j2,j1,-1
      if(aa(i2).eq.0)then
      aa(i2)=1
      bb(kk)=i2
      kk=kk-1
      lin=i2
      j2=lin-1
      goto 21
      endif
      enddo
      auxlin(1:srec)='gen21_9'
      call messag(1,0,0,0)
      endif
      goto 21
      endif
      if(lin.ne.n)then
      goto 17
      endif
      if(mflag(10).ne.0)then
      do i1=rhop1,n
      do i2=i1+1,n
      if(xg(i1,i2).ne.0)then
      if(aa(i1)+aa(i2).ne.0)then
      if(mflag(10).gt.0)then
      goto 05
      else
      goto 62
      endif
      endif
      endif
      enddo
      if(xg(i1,i1).ne.0)then
      if(mflag(10).gt.0)then
      goto 05
      else
      goto 62
      endif
      endif
      enddo
      if(mflag(10).lt.0)then
      goto 05
      endif
      endif
   62 continue
      if(mflag(9).lt.0)then
      do i1=rhop1,n
      do i2=i1,n
      if(xg(i1,i2).gt.1)then
      goto 34
      endif
      enddo
      enddo
      goto 05
      endif
   34 continue
      do i1=1,n
      do i2=1,i1-1
      j1=xg(i2,i1)
      g(i1,i2)=j1
      g(i2,i1)=j1
      enddo
      g(i1,i1)=xg(i1,i1)
      enddo
      if(mflag(1).ne.0)then
      call umpi(1,ii)
      if(ii.ne.mflag(1))then
      goto 05
      endif
      endif
      if(mflag(2).ne.0)then
      call umpi(2,ii)
      if(ii.ne.mflag(2))then
      goto 05
      endif
      endif
      if(mflag(4).ne.0)then
      call umpi(4,ii)
      if(ii.ne.mflag(4))then
      goto 05
      endif
      endif
      if(mflag(6).ne.0)then
      call umvi(2,ii)
      if(ii.ne.mflag(6))then
      goto 05
      endif
      endif
      if(mflag(8).ne.0)then
      if(rho(1).ne.1)then
      call umvi(1,ii)
      else
      call umvi(2,ii)
      endif
      if(ii.ne.mflag(8))then
      goto 05
      endif
      endif
      ntadp=0
      if(mflag(19).ne.0)then
      call umpi(3,ii)
      if(ii.ne.1)then
      auxlin(1:srec)='gen21_11'
      call messag(1,0,0,0)
      endif
      endif
      nsym=0
      do i1=1,rho(1)-1
      do i2=i1+1,rho(1)
      p1s(i1,i2)=0
      enddo
      enddo
      do i1=1,n
      xp(i1)=i1
      enddo
      goto 93
   77 continue
      do i1=xset(n),1,-1
      do i2=uset(i1)-1,uset(i1-1)+1,-1
      if(xp(i2).lt.xp(i2+1))then
      goto 102
      endif
      enddo
      enddo
      goto 63
  102 continue
      j1=uset(i1)
      do i1=i2+2,j1
      if(xp(i1).lt.xp(i2))then
      goto 202
      endif
      enddo
      i1=j1+1
  202 continue
      i1=i1-1
      ii=xp(i1)
      xp(i1)=xp(i2)
      xp(i2)=ii
      i1=i2+1
      i2=j1
  204 continue
      if(i1.lt.i2)then
      ii=xp(i1)
      xp(i1)=xp(i2)
      xp(i2)=ii
      i1=i1+1
      i2=i2-1
      goto 204
      endif
      do i1=j1+1,n
      xp(i1)=i1
      enddo
   93 continue
      if(rho(1).gt.0)then
      if(xp(rho(1)).ne.rho(1))then
      goto 63
      endif
      endif
      do i1=rhop1,n-1
      do i2=i1,n
      ii=g(xp(i1),xp(i2))-g(i1,i2)
      if(ii.gt.0)then
      goto 05
      endif
      if(ii.lt.0)then
      j1=xset(i2)
      j3=uset(j1)
      do i3=i2+1,j3-1
      do i4=i3+1,j3
      if(xp(i3).lt.xp(i4))then
      ii=xp(i3)
      xp(i3)=xp(i4)
      xp(i4)=ii
      endif
      enddo
      enddo
      do i3=j1+1,xset(n)
      j2=j3+1
      j3=uset(i3)
      j4=j2+j3
      do i4=j2,j3
      xp(i4)=j4-i4
      enddo
      enddo
      goto 77
      endif
      enddo
      enddo
      if(rho(1).eq.0)then
      goto 114
      endif
      i1=1
  110 continue
      j1=vmap(i1,1)
      if(xp(j1).eq.j1)then
      i1=i1+xn(j1)
      if(i1.le.rho(1))then
      goto 110
      endif
      goto 114
      else
      p1s(i1,a1(xp(j1)))=1
      goto 77
      endif
  114 continue
      if(jflag(11).eq.0)then
      if(psym(0).lt.0)then
      psyms=0
      psym(0)=stibs(1)
      call aoib(1)
      stib(stibs(1))=eoa
      endif
      jj=n-rho(1)
      ii=nsym*jj-psyms
      if(ii.gt.0)then
      if(psym(0).eq.stibs(1)-1-psyms)then
      call aoib(ii)
      psyms=psyms+ii
      stib(stibs(1))=eoa
      else
      auxlin(1:srec)='gen21_13'
      call messag(1,0,0,0)
      endif
      endif
      ii=nsym*jj
      if(ii.le.psyms)then
      if(nsym.gt.0)then
      ii=psym(0)+(ii-n)
      do i1=rhop1,n
      stib(ii+i1)=xp(i1)
      enddo
      endif
      else
      auxlin(1:srec)='gen21_14'
      call messag(1,0,0,0)
      endif
      endif
      nsym=nsym+1
      goto 77
   63 continue
      ns1=0
      do i1=2,rho(1)
      j1=i1-1
      if(vmap(j1,1).eq.vmap(i1,1))then
      p1s(j1,i1)=1
      endif
      do i2=1,j1
      if(p1s(i2,i1).eq.1)then
      ns1=ns1+1
      p1l(ns1)=i2
      p1r(ns1)=i1
      endif
      enddo
      enddo
      if(jflag(1)+jflag(10).eq.0)then
      goto 90
      endif
      do i1=1,rho(1)
      tail(i1)=i1
      head(i1)=vmap(i1,1)
      enddo
      do i1=1,degree(n)
      nemul(i1)=0
      enddo
      ii=rhop1
      do i1=rhop1,n
      do i2=i1,n
      jj=g(i1,i2)
      if(jj.eq.0)then
      eg(i1,i2)=0
      else
      eg(i1,i2)=ii
      if(i1.eq.i2)then
      jj=jj/2
      else
      nemul(jj)=nemul(jj)+1
      emul(jj,nemul(jj))=ii
      endif
      ii=ii+jj
      do i3=ii-jj,ii-1
      tail(i3)=i1
      head(i3)=i2
      enddo
      endif
      enddo
      enddo
      if(jflag(10).eq.0)then
      goto 90
      endif
      do i1=rhop1,nli
      intree(i1)=0
      enddo
      do i1=1,n
      bb(i1)=0
      enddo
      do i1=rhop1,n
      aa(i1)=0
      enddo
      kk=0
      ntree=0
      ii=degree(n)
   51 continue
      if(ii.le.0)then
      goto 52
      endif
      jj=nemul(ii)
   58 continue
      if(jj.le.0)then
      ii=ii-1
      goto 51
      endif
      j1=tail(emul(ii,jj))
      j2=head(emul(ii,jj))
      if(aa(j1).eq.aa(j2))then
      if(aa(j1).eq.0)then
      ntree=ntree+1
      aa(j1)=ntree
      aa(j2)=ntree
      else
      jj=jj-1
      goto 58
      endif
      else
      if(aa(j1).eq.0)then
      aa(j1)=aa(j2)
      elseif(aa(j2).eq.0)then
      aa(j2)=aa(j1)
      else
      j3=aa(j2)
      do i1=rhop1,n
      if(aa(i1).eq.j3)then
      aa(i1)=aa(j1)
      endif
      enddo
      endif
      endif
      intree(eg(j1,j2))=1
      bb(j1)=bb(j1)+1
      bb(j2)=bb(j2)+1
      kk=kk+1
      if(kk.lt.n-rhop1)then
      jj=jj-1
      goto 58
      endif
   52 continue
      do i1=1,rho(1)
      intree(i1)=1
      enddo
      do i1=1,nli
      do i2=1,rho(1)+loop
      flow(i1,i2)=0
      enddo
      enddo
      do i1=1,rho(1)
      flow(i1,i1)=1
      enddo
      ii=rhop1
      do i1=rhop1,nli
      if(intree(i1).eq.0)then
      flow(i1,ii)=1
      ii=ii+1
      endif
      enddo
      ii=rhop1
   73 continue
      if(bb(ii).eq.1)then
      goto 71
      elseif(ii.lt.n)then
      ii=ii+1
      goto 73
      else
      goto 83
      endif
   71 continue
      do i1=rhop1,ii-1
      if(xg(i1,ii).ne.0)then
      j1=eg(i1,ii)
      if(intree(j1).ne.0)then
      goto 75
      endif
      endif
      enddo
      do i1=ii+1,n
      if(xg(ii,i1).ne.0)then
      j1=eg(ii,i1)
      if(intree(j1).ne.0)then
      goto 75
      endif
      endif
      enddo
      auxlin(1:srec)='gen21_16'
      call messag(1,0,0,0)
   75 continue
      if(rho(1).gt.1)then
      if(xn(ii).gt.0)then
      do i1=1,rho(1)
      if(vmap(i1,1).eq.ii)then
      if(tail(j1).eq.ii)then
      flow(j1,i1)=flow(j1,i1)+1
      else
      flow(j1,i1)=flow(j1,i1)-1
      endif
      endif
      enddo
      endif
      endif
      do i1=rhop1,n
      if(i1.ne.ii)then
      if(i1.lt.ii)then
      jj=eg(i1,ii)
      else
      jj=eg(ii,i1)
      endif
      do i2=jj,jj+g(ii,i1)-1
      if(intree(i2).eq.0)then
      if((head(j1).eq.head(i2)).or.(tail(j1).eq.tail(i2)))then
      do i3=1,rho(1)+loop
      flow(j1,i3)=flow(j1,i3)-flow(i2,i3)
      enddo
      else
      do i3=1,rho(1)+loop
      flow(j1,i3)=flow(j1,i3)+flow(i2,i3)
      enddo
      endif
      endif
      enddo
      endif
      enddo
      intree(j1)=0
      bb(tail(j1))=bb(tail(j1))-1
      bb(head(j1))=bb(head(j1))-1
      ii=rhop1
      goto 73
   83 continue
      do i1=rhop1,nli
      ii=0
      jj=0
      do i2=1,rho(1)
      if(flow(i1,i2).eq.1)then
      ii=ii+1
      elseif(flow(i1,i2).eq.-1)then
      jj=jj+1
      endif
      enddo
      if((ii.gt.0).and.(jj.gt.0))then
      auxlin(1:srec)='gen21_17'
      call messag(1,0,0,0)
      endif
      if(2*(ii+jj).gt.rho(1))then
      if(ii.gt.jj)then
      do i2=1,rho(1)
      flow(i1,i2)=flow(i1,i2)-1
      enddo
      else
      do i2=1,rho(1)
      flow(i1,i2)=flow(i1,i2)+1
      enddo
      endif
      endif
      enddo
      nbri=0
      nrbri=0
      nsbri=0
      do i1=rhop1,nli
      do i2=rhop1,rho(1)+loop
      if(flow(i1,i2).ne.0)then
      flow(i1,0)=1
      goto 60
      endif
      enddo
      nbri=nbri+1
      do i2=1,rho(1)
      if(flow(i1,i2).ne.0)then
      if(mflag(3).gt.0)then
      goto 05
      endif
      flow(i1,0)=2
      nrbri=nrbri+1
      goto 60
      endif
      enddo
      nsbri=nsbri+1
      flow(i1,0)=3
   60 continue
      enddo
      if(mflag(3).lt.0)then
      if(nrbri.eq.0)then
      goto 05
      endif
      endif
      if(zbri(nbri).eq.0)then
      goto 05
      elseif(rbri(nrbri).eq.0)then
      goto 05
      elseif(sbri(nsbri).eq.0)then
      goto 05
      endif
      if(zcho(nli-rho(1)-nbri).eq.0)then
      goto 05
      endif
      if(mflag(5).ne.0)then
      do i1=rhop1,nli
      do i2=1,i1-1
      do i3=0,1
      do i4=rhop1,rho(1)+loop
      if(i3.eq.0)then
      ii=flow(i1,i4)+flow(i2,i4)
      else
      ii=flow(i1,i4)-flow(i2,i4)
      endif
      if(ii.ne.0)then
      goto 105
      endif
      enddo
      do i4=2,rho(1)
      ii=flow(i1,i4-1)-flow(i1,i4)
      jj=flow(i2,i4-1)-flow(i2,i4)
      if(i3.eq.0)then
      ii=ii+jj
      else
      ii=ii-jj
      endif
      if(ii.ne.0)then
      goto 105
      endif
      enddo
      if(mflag(5).gt.0)then
      goto 05
      else
      goto 74
      endif
 105  continue
      enddo
      enddo
      enddo
      if(mflag(5).lt.0)then
      goto 05
      endif
      endif
   74 continue
      do i1=1,nli
      do i2=1,rho(1)+loop
      if(abs(flow(i1,i2)).gt.1)then
      auxlin(1:srec)='gen21_18'
      call messag(1,0,0,0)
      endif
      enddo
      enddo
      do i1=1,rho(1)
      eg(i1,vmap(i1,1))=i1
      enddo
      do i1=rhop1,n
      ii=0
      do i2=1,rho(1)
      jj=0
      do i3=1,n
      if(i3.ne.i1)then
      j2=g(i1,i3)
      if(j2.gt.0)then
      if(i1.lt.i3)then
      j1=eg(i1,i3)
      else
      j1=eg(i3,i1)
      endif
      j3=head(j1)-i1
      do i4=j1,j1+j2-1
      if(j3.eq.0)then
      jj=jj+flow(i4,i2)
      else
      jj=jj-flow(i4,i2)
      endif
      enddo
      endif
      endif
      enddo
      if(jj.ne.ii)then
      if(i2.eq.1)then
      ii=jj
      else
      auxlin(1:srec)='gen21_19'
      call messag(1,0,0,0)
      endif
      endif
      enddo
      enddo
      do i1=rhop1,n
      do i2=rhop1,rho(1)+loop
      ii=0
      do i3=1,n
      if(i3.ne.i1)then
      j2=g(i1,i3)
      if(j2.gt.0)then
      if(i1.lt.i3)then
      j1=eg(i1,i3)
      else
      j1=eg(i3,i1)
      endif
      j3=head(j1)-i1
      do i4=j1,j1+j2-1
      if(j3.eq.0)then
      ii=ii+flow(i4,i2)
      else
      ii=ii-flow(i4,i2)
      endif
      enddo
      endif
      endif
      enddo
      if(ii.ne.0)then
      auxlin(1:srec)='gen21_20'
      call messag(1,0,0,0)
      endif
      enddo
      enddo
      goto 90
   97 continue
      auxlin(1:srec)='gen21_21'
      call messag(1,0,0,0)
   80 continue
      cntr21=0
   90 continue
      return
      end
      subroutine init0
      implicit integer(a-z)
      save
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( sxbuff=2040 )
      parameter ( srec=81, ssrec=62 )
      common/z12g/jflag(1:11),mflag(1:24)
      character*(srec) auxlin
      common/z22g/auxlin
      character*(sxbuff) stxb
      common/z27g/stxb
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z48g/iref,wiref
      common/z51g/aplus,aminus,azero,anine,aeq
      common/z52g/acf0(0:127),acf1(0:127),lfeed,aspace,dquote,squote
      wiref=4
      iref=1
      do i1=1,wiref
      iref=10*iref
      enddo
      iref=iref-1
      do i1=1,11
      jflag(i1)=0
      enddo
      do i1=1,24
      mflag(i1)=0
      enddo
      aplus=ichar('+')
      aminus=ichar('-')
      azero=ichar('0')
      anine=ichar('9')
      aeq=ichar('=')
      lfeed=10
      aspace=ichar(' ')
      dquote=ichar('"')
      squote=ichar("'")
      do i1=0,31
      acf0(i1)=-1
      enddo
      do i1=32,126
      acf0(i1)=0
      enddo
      acf0(127)=-1
      acf0(lfeed)=1
      acf0(aspace)=2
      acf0(dquote)=3
      acf0(squote)=4
      acf0(ichar('('))=5
      acf0(ichar(')'))=6
      acf0(ichar(','))=7
      acf0(ichar(';'))=8
      acf0(aeq)=9
      acf0(ichar('['))=10
      acf0(ichar(']'))=11
      do i1=0,127
      acf1(i1)=-1
      enddo
      do i1=azero,anine
      acf1(i1)=1
      enddo
      do i1=65,90
      acf1(i1)=2
      enddo
      acf1(95)=0
      do i1=97,122
      acf1(i1)=3
      enddo
      if(sxbuff.lt.2040)then
      call uput(5)
      endif
      if(sibuff.lt.8184)then
      call uput(6)
      endif
      if(scbuff.lt.8184)then
      call uput(7)
      endif
      if((mimal.lt.1).or.(mimal.gt.2))then
      call uput(8)
      endif
      if((mcmal.lt.1).or.(mcmal.gt.2))then
      call uput(9)
      endif
      if(srec.lt.81)then
      call uput(10)
      endif
      if((ssrec.ge.min(srec,128)).or.(ssrec.le.48))then
      call uput(11)
      endif
      imal(0)=1
      do i1=0,mimal
      stibs(i1)=0
      enddo
      cmal(0)=1
      do i1=0,mcmal
      stcbs(i1)=0
      enddo
      call init1
      return
      end
      subroutine init1
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( sxbuff=2040 )
      parameter ( eoa=-127, nap=-2047 )
      common/z5in/dunit,munit,ounit,sunit,funit
      common/z5g/psym(0:0),psyms,nsym
      common/z13g/kloo(15:32,1:4),loopt(11:19)
      common/z21g/punct1(0:0)
      character*(srec) auxlin
      common/z22g/auxlin
      common/z26g/kes(0:0),kle(0:0),pstke(0:0),wstke(0:0)
      character*(sxbuff) stxb
      common/z27g/stxb
      common/z28g/wera(0:0),werb(0:0),nwer
      common/z29g/pkey(0:0),wkey(0:0),prevl(0:0),nextl(0:0),popt3(0:0),
     :wopt3(0:0),fopt3(0:0),vopt3(0:0),popt5(0:0),wopt5(0:0),copt5(0:0)
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z32g/sucpal(0:11,0:11)
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z44g/xtstrp(0:0),xtstrl(0:0)
      common/z45g/qdatp,qdatl,qvp,qvl,dprefp,dprefl
      common/z46g/popt1(0:0),wopt1(0:0),copt1(0:0),popt7(0:0),
     :wopt7(0:0),copt7(0:0)
      common/z47g/ndiagp,ndiagl,hhp,hhl,doffp,doffl,noffp,noffl,wsint
      common/z48g/iref,wiref
      common/z49g/popt9(0:0),wopt9(0:0),copt9(0:0)
      common/z50g/tfta(0:0),tftb(0:0),tftic(0:0),ntft
      common/z51g/aplus,aminus,azero,anine,aeq
      common/z52g/acf0(0:127),acf1(0:127),lfeed,aspace,dquote,squote
      integer weri(0:0)
      j1=0
      do i1=1,srec
      j1=j1+2
      stxb(j1-1:j1)='00'
      enddo
      j1=j1+1
      stxb(j1:j1)=';'
      call spak(stxb(1:j1),jj,0,0,0)
      wsint=15
      jj=wsint+1
      ndiagp=stcbs(1)
      call aocb(jj)
      ii=stcbs(1)
      stcb(ii:ii)=char(lfeed)
      call sdiag(3,-1)
      hhp=stcbs(1)
      call aocb(jj)
      ii=stcbs(1)
      stcb(ii:ii)=char(lfeed)
      noffp=stcbs(1)
      noffl=1
      call aocb(jj)
      ii=noffp+1
      stcb(ii:ii)=char(azero)
      ii=stcbs(1)
      stcb(ii:ii)=char(lfeed)
      doffp=stcbs(1)
      doffl=1
      call aocb(jj)
      ii=doffp+1
      stcb(ii:ii)=char(azero)
      ii=stcbs(1)
      stcb(ii:ii)=char(lfeed)
      ii=0
      call spak('817182657013191419;',ii,1,0,0)
      qvp=stib(stibs(1)+1)
      qvl=stib(stibs(1)+2)
      dunit=0
      munit=0
      ounit=0
      sunit=0
      funit=0
      call spak('817182657014686584;',ii,1,0,0)
      qdatp=stib(stibs(1)+1)
      qdatl=stib(stibs(1)+2)
      call spak('6885657613;',ii,1,0,0)
      dprefp=stib(stibs(1)+1)
      dprefl=stib(stibs(1)+2)
      nxts=0
      call spak('78797813808279806571658479828326000008;',nxts,0,0,0)
      call spak('808279806571658479828326000008;',nxts,0,0,0)
      call spak('4611;',nxts,0,0,0)
      call spak('4613;',nxts,0,0,0)
      call spak('3511;',nxts,0,0,0)
      call spak('3513;',nxts,0,0,0)
      call spak('866982847367698326000008;',nxts,0,0,0)
      call spak('861368697182696983;',nxts,0,0,0)
      call spak('036873657182657783;',nxts,0,0,0)
      call spak('03718265807283;',nxts,0,0,0)
      call spak('00000000000000847984657600290000;',nxts,0,0,0)
      call spak('006873657182657783;',nxts,0,0,0)
      call spak('00718265807283;',nxts,0,0,0)
      call spak('0900;',nxts,0,0,0)
      call spak('14141414;',nxts,0,0,0)
      call spak('0000001010;',nxts,0,0,0)
      call spak('000000876582787378712600;',nxts,0,0,0)
      call spak('00000069828279822600;',nxts,0,0,0)
      call spak('7073766900;',nxts,0,0,0)
      call spak('77796869761370737669;',nxts,0,0,0)
      call spak('767366826582891370737669;',nxts,0,0,0)
      call spak('83848976691370737669;',nxts,0,0,0)
      call spak('7985848085841370737669;',nxts,0,0,0)
      call spak('12007673786900;',nxts,0,0,0)
      call spak('1200767378698300;',nxts,0,0,0)
      call spak('101010;',nxts,0,0,0)
      kk=2
      call aoib(kk)
      jj=stibs(1)
      do i1=1,kk
      stib(jj)=eoa
      jj=jj-1
      enddo
      ii=nxts+1
      call trm(kk,ii)
      xtstrl(0)=stibs(1)-ii
      xtstrp(0)=xtstrl(0)-ii
      psym(0)=-1
      nos=0
      uc=1
      nkey=0
      call spak('798584808584,0,1;',nkey,0,uc,nos)
      call spak('8384897669,1,2;',nkey,0,uc,nos)
      call spak('7779686976,2,4;',nkey,0,uc,nos)
      call spak('7378,4,5;',nkey,0,uc,nos)
      call spak('798584,5,6;',nkey,0,uc,nos)
      call spak('7679798083,6,7;',nkey,0,uc,nos)
      call spak('76797980637779776978848577,7,8;',nkey,0,uc,nos)
      call spak('79808473797883,8,9;',nkey,0,uc,nos)
      call spak('737868698863797070836984,9,-1;',nkey,0,uc,nos)
      kk=4
      call aoib(kk)
      jj=stibs(1)
      do i1=1,kk
      stib(jj)=eoa
      jj=jj-1
      enddo
      ii=nkey+1
      call trm(kk,ii)
      nextl(0)=stibs(1)-ii
      prevl(0)=nextl(0)-ii
      wkey(0)=prevl(0)-ii
      pkey(0)=wkey(0)-ii
      nos=0
      uc=1
      nopt1=0
      call spak('84828569,1;',nopt1,0,uc,nos)
      call spak('7065768369,-1;',nopt1,0,uc,nos)
      kk=3
      call aoib(kk)
      jj=stibs(1)
      do i1=1,kk
      stib(jj)=eoa
      jj=jj-1
      enddo
      ii=nopt1+1
      call trm(kk,ii)
      copt1(0)=stibs(1)-ii
      wopt1(0)=copt1(0)-ii
      popt1(0)=wopt1(0)-ii
      nos=1
      uc=1
      nopt3=0
      call spak('787966827368716983,1,1;',nopt3,0,uc,nos)
      call spak('7978698073,1,1;',nopt3,0,uc,nos)
      call spak('66827368716983,1,-1;',nopt3,0,uc,nos)
      call spak('7978698082,1,-1;',nopt3,0,uc,nos)
      call spak('78798366827368716983,2,1;',nopt3,0,uc,nos)
      call spak('78798465688079766983,2,1;',nopt3,0,uc,nos)
      call spak('8366827368716983,2,-1;',nopt3,0,uc,nos)
      call spak('8465688079766983,2,-1;',nopt3,0,uc,nos)
      call spak('78798266827368716983,3,1;',nopt3,0,uc,nos)
      call spak('8266827368716983,3,-1;',nopt3,0,uc,nos)
      call spak('79788372697676,4,1;',nopt3,0,uc,nos)
      call spak('7970708372697676,4,-1;',nopt3,0,uc,nos)
      call spak('7879837371776583,5,1;',nopt3,0,uc,nos)
      call spak('837371776583,5,-1;',nopt3,0,uc,nos)
      call spak('7978837269767688,6,1;',nopt3,0,uc,nos)
      call spak('797070837269767688,6,-1;',nopt3,0,uc,nos)
      call spak('6789677673,7,1;',nopt3,0,uc,nos)
      call spak('7879786789677673,7,-1;',nopt3,0,uc,nos)
      call spak('7879837865737683,8,1;',nopt3,0,uc,nos)
      call spak('837865737683,8,-1;',nopt3,0,uc,nos)
      call spak('837377807669,9,1;',nopt3,0,uc,nos)
      call spak('787984837377807669,9,-1;',nopt3,0,uc,nos)
      call spak('667380658284,10,1;',nopt3,0,uc,nos)
      call spak('787978667380658284,10,-1;',nopt3,0,uc,nos)
      call spak('8479807976,11,1;',nopt3,0,uc,nos)
      call spak('707679798083,12,1;',nopt3,0,uc,nos)
      call spak('787984707679798083,12,-1;',nopt3,0,uc,nos)
      call spak('7879798084,13,1;',nopt3,0,uc,nos)
      kk=4
      call aoib(kk)
      jj=stibs(1)
      do i1=1,kk
      stib(jj)=eoa
      jj=jj-1
      enddo
      ii=nopt3+1
      call trm(kk,ii)
      vopt3(0)=stibs(1)-ii
      fopt3(0)=vopt3(0)-ii
      wopt3(0)=fopt3(0)-ii
      popt3(0)=wopt3(0)-ii
      nos=0
      uc=1
      nopt5=0
      call spak('7380827980,1;',nopt5,0,uc,nos)
      call spak('668273687169,2;',nopt5,0,uc,nos)
      call spak('6772798268,3;',nopt5,0,uc,nos)
      call spak('82668273687169,4;',nopt5,0,uc,nos)
      call spak('83668273687169,5;',nopt5,0,uc,nos)
      call spak('86838577,6;',nopt5,0,uc,nos)
      call spak('80838577,7;',nopt5,0,uc,nos)
      call spak('6976737875,11;',nopt5,0,uc,nos)
      call spak('8076737875,12;',nopt5,0,uc,nos)
      call spak('8676737875,-1;',nopt5,0,uc,nos)
      kk=3
      call aoib(kk)
      jj=stibs(1)
      do i1=1,kk
      stib(jj)=eoa
      jj=jj-1
      enddo
      ii=nopt5+1
      call trm(kk,ii)
      copt5(0)=stibs(1)-ii
      wopt5(0)=copt5(0)-ii
      popt5(0)=wopt5(0)-ii
      kk=0
      do i1=1,nopt5
      j1=stib(copt5(0)+i1)
      if(j1.ne.-1)then
      if(j1.gt.kk)then
      ntft=ntft+1
      kk=j1
      elseif(j1.lt.kk)then
      auxlin(1:srec)='init1_1'
      call messag(1,0,0,0)
      endif
      endif
      enddo
      ii=ntft+1
      jj=2*ii
      call aoib(jj)
      tftb(0)=stibs(1)-ii
      tfta(0)=tftb(0)-ii
      stib(tftb(0))=eoa
      stib(stibs(1))=eoa
      jj=kk+1
      tftic(0)=stibs(1)
      call aoib(jj)
      stib(stibs(1))=eoa
      do i1=tftic(0)+1,tftic(0)+kk
      stib(tftic(0)+i1)=0
      enddo
      ii=0
      do i1=1,nopt5
      j1=tftic(0)+stib(copt5(0)+i1)
      if(stib(j1).eq.0)then
      ii=ii+1
      stib(j1)=ii
      endif
      enddo
      nopt7=0
      nos=1
      uc=1
      call spak('6988846982786576,5;',nopt7,0,uc,nos)
      call spak('78798465688079766983,1;',nopt7,0,uc,nos)
      kk=3
      call aoib(kk)
      jj=stibs(1)
      do i1=1,kk
      stib(jj)=eoa
      jj=jj-1
      enddo
      ii=nopt7+1
      call trm(kk,ii)
      copt7(0)=stibs(1)-ii
      wopt7(0)=copt7(0)-ii
      popt7(0)=wopt7(0)-ii
      nopt9=0
      nos=1
      uc=1
      call spak('69886776,1;',nopt9,0,uc,nos)
      call spak('73786776,0;',nopt9,0,uc,nos)
      kk=3
      call aoib(kk)
      jj=stibs(1)
      do i1=1,kk
      stib(jj)=eoa
      jj=jj-1
      enddo
      ii=nopt9+1
      call trm(kk,ii)
      copt9(0)=stibs(1)-ii
      wopt9(0)=copt9(0)-ii
      popt9(0)=wopt9(0)-ii
      nos=0
      uc=0
      nstke=0
      call spak('8082797679718569,1,0;',nstke,0,uc,nos)
      call spak('68736571826577,2,0;',nstke,0,uc,nos)
      call spak('6980737679718569,3,0;',nstke,0,uc,nos)
      call spak('69887384,4,0;',nstke,0,uc,nos)
      call spak('677977776578686376797980,11,5;',nstke,0,uc,nos)
      call spak('6779777765786863767378696376797980,12,5;',nstke,0,
     :uc,nos)
      call spak('73786376797980,13,2;',nstke,0,uc,nos)
      call spak('7985846376797980,14,2;',nstke,0,uc,nos)
      call spak('808279806571658479826376797980,15,2;',nstke,0,uc,nos)
      call spak('8669828469886376797980,16,2;',nstke,0,uc,nos)
      call spak('8265896376797980,17,2;',nstke,0,uc,nos)
      call spak('697868,19,7;',nstke,0,uc,nos)
      call spak('66656775,21,7;',nstke,0,uc,nos)
      call spak('78698776737869,22,7;',nstke,0,uc,nos)
      call spak('7073697668,31,2;',nstke,0,uc,nos)
      kloo(15,1)=1
      kloo(15,2)=1
      kloo(15,3)=0
      kloo(15,4)=1
      call spak('7779776978848577,32,2;',nstke,0,uc,nos)
      kloo(16,1)=1
      kloo(16,2)=1
      kloo(16,3)=0
      kloo(16,4)=1
      call spak('68856576137073697668,33,2;',nstke,0,uc,nos)
      kloo(17,1)=1
      kloo(17,2)=1
      kloo(17,3)=0
      kloo(17,4)=1
      call spak('68856576137779776978848577,34,2;',nstke,0,uc,nos)
      kloo(18,1)=1
      kloo(18,2)=1
      kloo(18,3)=0
      kloo(18,4)=1
      call spak('70736976686383737178,35,2;',nstke,0,uc,nos)
      kloo(19,1)=1
      kloo(19,2)=1
      kloo(19,3)=0
      kloo(19,4)=1
      call spak('70736976686384898069,40,2;',nstke,0,uc,nos)
      kloo(20,1)=1
      kloo(20,2)=1
      kloo(20,3)=0
      kloo(20,4)=1
      call spak('80827980657165847982637378686988,41,2;',nstke,0,uc,nos)
      kloo(21,1)=0
      kloo(21,2)=1
      kloo(21,3)=0
      kloo(21,4)=1
      call spak('7073697668637378686988,42,2;',nstke,0,uc,nos)
      kloo(22,1)=1
      kloo(22,2)=1
      kloo(22,3)=0
      kloo(22,4)=1
      call spak('826589637378686988,43,2;',nstke,0,uc,nos)
      kloo(23,1)=1
      kloo(23,2)=1
      kloo(23,3)=0
      kloo(23,4)=1
      call spak('866982846988637378686988,44,2;',nstke,0,uc,nos)
      kloo(24,1)=1
      kloo(24,2)=1
      kloo(24,3)=1
      kloo(24,4)=1
      call spak('68856576137073697668637378686988,45,2;',nstke,0,uc,nos)
      kloo(25,1)=0
      kloo(25,2)=1
      kloo(25,3)=0
      kloo(25,4)=1
      call spak('6885657613826589637378686988,46,2;',nstke,0,uc,nos)
      kloo(26,1)=0
      kloo(26,2)=1
      kloo(26,3)=0
      kloo(26,4)=1
      call spak('6885657613866982846988637378686988,47,2;',nstke,0,
     :uc,nos)
      kloo(27,1)=0
      kloo(27,2)=1
      kloo(27,3)=0
      kloo(27,4)=1
      call spak('86698284698863686971826969,48,2;',nstke,0,uc,nos)
      kloo(28,1)=1
      kloo(28,2)=1
      kloo(28,3)=1
      kloo(28,4)=1
      call spak('688565761386698284698863686971826969,49,2;',nstke,0,
     :uc,nos)
      kloo(29,1)=0
      kloo(29,2)=1
      kloo(29,3)=0
      kloo(29,4)=1
      call spak('766971637378686988,51,2;',nstke,0,uc,nos)
      kloo(30,1)=1
      kloo(30,2)=0
      kloo(30,3)=0
      kloo(30,4)=0
      call spak('7378637378686988,52,2;',nstke,0,uc,nos)
      kloo(31,1)=1
      kloo(31,2)=0
      kloo(31,3)=0
      kloo(31,4)=0
      call spak('798584637378686988,53,2;',nstke,0,uc,nos)
      kloo(32,1)=1
      kloo(32,2)=0
      kloo(32,3)=0
      kloo(32,4)=0
      call spak('83737178,61,2;',nstke,0,uc,nos)
      call spak('7773788583,62,2;',nstke,0,uc,nos)
      call spak('68736571826577637378686988,71,6;',nstke,0,uc,nos)
      call spak('838977776984828963706567847982,72,2;',nstke,0,uc,nos)
      call spak('838977776984828963788577666982,73,2;',nstke,0,uc,nos)
      call spak('8082798065716584798283,74,2;',nstke,0,uc,nos)
      call spak('76697183,75,2;',nstke,0,uc,nos)
      call spak('7679798083,76,2;',nstke,0,uc,nos)
      call spak('8669828473676983,77,2;',nstke,0,uc,nos)
      call spak('76697183637378,78,2;',nstke,0,uc,nos)
      call spak('7669718363798584,79,2;',nstke,0,uc,nos)
      call spak('80827971826577,81,5;',nstke,0,uc,nos)
      call spak('677977776578686368658465,82,5;',nstke,0,uc,nos)
      kk=4
      call aoib(kk)
      jj=stibs(1)
      do i1=1,kk
      stib(jj)=eoa
      jj=jj-1
      enddo
      ii=nstke+1
      call trm(kk,ii)
      kle(0)=stibs(1)-ii
      kes(0)=kle(0)-ii
      wstke(0)=kes(0)-ii
      pstke(0)=wstke(0)-ii
      loopt(11)=1
      loopt(12)=2
      loopt(13)=3
      loopt(14)=3
      loopt(15)=4
      loopt(16)=5
      loopt(17)=6
      loopt(18)=0
      loopt(19)=-1
      punct1(0)=stibs(1)+1
      ii=128
      call aoib(ii)
      do i1=0,127
      stib(punct1(0)+i1)=0
      enddo
      stib(punct1(0)+ichar('='))=1
      stib(punct1(0)+ichar(','))=2
      stib(punct1(0)+ichar(';'))=3
      stib(punct1(0)+squote)=4
      stib(punct1(0)+ichar('['))=5
      stib(punct1(0)+ichar(']'))=6
      stib(punct1(0)+ichar('('))=-7
      stib(punct1(0)+ichar(')'))=7
      do i1=0,11
      do i2=0,11
      sucpal(i1,i2)=-1
      enddo
      enddo
      sucpal(0,2)=0
      sucpal(0,10)=1
      sucpal(1,0)=1
      sucpal(1,1)=1
      sucpal(1,2)=1
      sucpal(1,7)=1
      sucpal(1,8)=2
      sucpal(1,9)=3
      sucpal(1,11)=10
      sucpal(2,0)=2
      sucpal(2,1)=2
      sucpal(2,2)=2
      sucpal(2,9)=3
      sucpal(3,0)=3
      sucpal(3,1)=3
      sucpal(3,2)=3
      sucpal(3,3)=-1
      sucpal(3,4)=5
      sucpal(3,5)=6
      sucpal(3,7)=2
      sucpal(3,11)=10
      sucpal(4,0)=4
      sucpal(4,2)=4
      sucpal(4,3)=3
      sucpal(4,4)=4
      sucpal(4,5)=4
      sucpal(4,6)=4
      sucpal(4,7)=4
      sucpal(4,8)=4
      sucpal(4,9)=4
      sucpal(4,10)=4
      sucpal(4,11)=4
      sucpal(5,0)=5
      sucpal(5,2)=5
      sucpal(5,3)=5
      sucpal(5,4)=3
      sucpal(5,5)=5
      sucpal(5,6)=5
      sucpal(5,7)=5
      sucpal(5,8)=5
      sucpal(5,9)=5
      sucpal(5,10)=5
      sucpal(5,11)=5
      sucpal(6,0)=6
      sucpal(6,1)=6
      sucpal(6,2)=6
      sucpal(6,3)=-1
      sucpal(6,4)=8
      sucpal(6,6)=9
      sucpal(6,7)=6
      sucpal(7,0)=7
      sucpal(7,2)=7
      sucpal(7,3)=6
      sucpal(7,4)=7
      sucpal(7,5)=7
      sucpal(7,6)=7
      sucpal(7,7)=7
      sucpal(7,8)=7
      sucpal(7,9)=7
      sucpal(7,10)=7
      sucpal(7,11)=7
      sucpal(8,0)=8
      sucpal(8,2)=8
      sucpal(8,3)=8
      sucpal(8,4)=6
      sucpal(8,5)=8
      sucpal(8,6)=8
      sucpal(8,7)=8
      sucpal(8,8)=8
      sucpal(8,9)=8
      sucpal(8,10)=8
      sucpal(8,11)=8
      sucpal(9,1)=9
      sucpal(9,2)=9
      sucpal(9,7)=2
      sucpal(9,11)=10
      sucpal(10,1)=11
      nwer=0
      call spak('7779686976006765780066690083807673840073788479006873'
     ://'8374797378840083856613777968697683,1;',nwer,0,0,0)
      call spak('6584007669658384007978690070736976680073830078798400'
     ://'806582840079700065788900866982846988,2;',nwer,0,0,0)
      call spak('7779686976137073766900726583006885807673676584690086'
     ://'69828473676983,3;',nwer,0,0,0)
      call spak('6977808489138384827378710067797883846578840068698469'
     ://'678469682600,4;',nwer,0,0,0)
      call spak('6977808489138384827378710070857867847379780068698469'
     ://'678469682600,5;',nwer,0,0,0)
      call spak('788576760070857867847379780068698469678469682600,6;',
     :nwer,0,0,0)
      call spak('7879008669828469880077658900767378751385800087738472'
     ://'00737867797773787100806582847367766900,7;',nwer,0,0,0)
      call spak('7879008669828469880077658900767378751385800087738472'
     ://'00798584717973787100806582847367766900,8;',nwer,0,0,0)
      call spak('6988846982786576008065828473677669830067657878798400'
     ://'666900677978786967846968,9;',nwer,0,0,0)
      call spak('7968680078857766698200797000698884698278657600657884'
     ://'7367797777858473787100707369766883,10;',nwer,0,0,0)
      call spak('67726967750079808473797883,11;',nwer,0,0,0)
      call spak('8482738673657600697673787500838465846977697884,12;',
     :nwer,0,0,0)
      call spak('6779788482656873678479828900697673787500838465846977'
     ://'697884,13;',nwer,0,0,0)
      call spak('0280767378750200826981857382698300687365718265778300'
     ://'8773847200658400766965838400180076697183,14;',nwer,0,0,0)
      kk=3
      call aoib(kk)
      jj=stibs(1)
      do i1=1,kk
      stib(jj)=eoa
      jj=jj-1
      enddo
      ii=nwer+1
      call trm(kk,ii)
      weri(0)=stibs(1)-ii
      werb(0)=weri(0)-ii
      wera(0)=werb(0)-ii
      do i1=1,nwer
      stib(werb(0)+i1)=stib(werb(0)+i1)-1+stib(wera(0)+i1)
      if(stib(weri(0)+i1).ne.i1)then
      auxlin(1:srec)='init1_2'
      call messag(1,0,0,0)
      endif
      enddo
      call aoib(-ii)
      weri(0)=0
      wiref=8
      iref=1
      do i1=1,wiref
      ii=iref
      iref=10*ii
      if(iref/ii.ne.10)then
      call uput(3)
      endif
      enddo
      iref=iref-1
      return
      end
      integer function lupty(lcode)
      implicit integer(a-z)
      save
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( srec=81, ssrec=62 )
      parameter ( loopt1=11, loopt2=19 )
      common/z13g/kloo(15:32,1:4),loopt(11:19)
      character*(srec) auxlin
      common/z22g/auxlin
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      lupty=0
      if((lcode.lt.loopt1).or.(lcode.gt.loopt2))then
      lupty=-1
      elseif(loopt(lcode).eq.0)then
      lupty=-1
      endif
      if(lupty.lt.0)then
      auxlin(1:srec)='lupty_1'
      call messag(1,0,0,0)
      endif
      lupty=loopt(lcode)
      return
      end
      subroutine aocb(delta)
      implicit integer(a-z)
      save
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( srec=81, ssrec=62 )
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      character*(srec) auxlin
      common/z22g/auxlin
      integer delta
      if(delta.gt.0)then
      if(delta.gt.scbuff.or.stcbs(1).gt.scbuff-delta)then
      call uput(2)
      endif
      elseif(delta.lt.0)then
      if(delta.lt.-stcbs(1))then
      auxlin(1:srec)='aocb_1'
      call messag(1,0,0,0)
      endif
      endif
      stcbs(1)=stcbs(1)+delta
      return
      end
      subroutine vaocb(delta)
      implicit integer(a-z)
      save
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( srec=81, ssrec=62 )
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      character*(srec) auxlin
      common/z22g/auxlin
      integer delta
      if((delta.lt.0).or.(stcbs(1).lt.0))then
      auxlin(1:srec)='vaocb_1'
      call messag(1,0,0,0)
      elseif(stcbs(1).gt.scbuff-delta)then
      call uput(2)
      endif
      return
      end
      subroutine aoib(delta)
      implicit integer(a-z)
      save
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( srec=81, ssrec=62 )
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      character*(srec) auxlin
      common/z22g/auxlin
      integer delta
      if(delta.gt.0)then
      if(delta.gt.sibuff.or.stibs(1).gt.sibuff-delta)then
      call uput(1)
      endif
      elseif(delta.lt.0)then
      if(delta.lt.-stibs(1))then
      auxlin(1:srec)='aoib_1'
      call messag(1,0,0,0)
      endif
      endif
      stibs(1)=stibs(1)+delta
      return
      end
      subroutine vaoib(delta)
      implicit integer(a-z)
      save
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( srec=81, ssrec=62 )
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      character*(srec) auxlin
      common/z22g/auxlin
      integer delta
      if((delta.lt.0).or.(stibs(1).lt.0))then
      auxlin(1:srec)='vaoib_1'
      call messag(1,0,0,0)
      elseif(stibs(1).gt.sibuff-delta)then
      call uput(1)
      endif
      return
      end
      subroutine mput(istop,nl1,nl2,nf1)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( sxbuff=2040 )
      common/z12g/jflag(1:11),mflag(1:24)
      character*(srec) auxlin
      common/z22g/auxlin
      character*(sxbuff) stxb
      common/z27g/stxb
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      common/z44g/xtstrp(0:0),xtstrl(0:0)
      common/z45g/qdatp,qdatl,qvp,qvl,dprefp,dprefl
      common/z48g/iref,wiref
      common/z51g/aplus,aminus,azero,anine,aeq
      common/z52g/acf0(0:127),acf1(0:127),lfeed,aspace,dquote,squote
      if(jflag(7).eq.0)then
      print *,''
      endif
      call hrul(2)
      print *,stxb(1:ssrec-1)
      slin1=ssrec-1
      slin2=2*srec
      slin3=(scbuff-wiref)-stcbs(1)
      stab=0
      do j2=srec-1,1,-1
      if(ichar(auxlin(j2:j2)).ne.aspace)then
      goto 30
      endif
      enddo
      j1=0
      j2=0
      goto 50
   30 continue
      do j1=1,j2
      if(ichar(auxlin(j1:j1)).ne.aspace)then
      goto 50
      endif
      enddo
   50 continue
      if(jflag(6).eq.0)then
      print *,'   error: '//auxlin(j1:j2)
      goto 70
      endif
      if(istop.eq.0)then
      j3=17
      else
      j3=18
      endif
      i1=stib(xtstrp(0)+j3)
      i2=stib(xtstrl(0)+j3)
      stab=i2-1
      if(i2.ge.slin2)then
      print *,''
      print *,' error: mput_1'
      print *,''
      stop
      endif
      stxb(1:i2)=stcb(i1:i1-1+i2)
      if(nf1.ge.0)then
      if(j1.gt.0)then
      j3=min(slin2-i2,j2-j1+1)
      if(j3.gt.0)then
      stxb(i2+1:i2+j3)=auxlin(j1:j1-1+j3)
      i2=i2+j3
      endif
      endif
      endif
      if(nf1.ne.0)then
      j4=abs(nf1)
      if(nf1.gt.0)then
      if(i2.lt.slin2)then
      i2=i2+1
      stxb(i2:i2)=stcb(1:1)
      endif
      endif
      j3=18+j4
      i=stib(xtstrl(0)+j3)
      if(i2+i.le.slin2)then
      j3=stib(xtstrp(0)+j3)
      stxb(i2+1:i2+i)=stcb(j3:j3-1+i)
      endif
      i2=i2+i
      if(j4.eq.1)then
      if(i2+qdatl.le.slin2)then
      stxb(i2+1:i2+qdatl)=stcb(qdatp:qdatp-1+qdatl)
      endif
      i2=i2+qdatl
      endif
      endif
      if(nf1.lt.0)then
      if(j1.gt.0)then
      if(j2-j1+2.le.slin2-i2)then
      stxb(i2+1:i2+2+(j2-j1))=stcb(1:1)//auxlin(j1:j2)
      i2=i2+(j2-j1)+2
      endif
      endif
      endif
      if(nl1.gt.0)then
      if(nl1.eq.nl2)then
      i=stib(xtstrl(0)+24)
      if(i2+i.le.slin2)then
      j3=stib(xtstrp(0)+24)
      stxb(i2+1:i2+i)=stcb(j3:j3-1+i)
      endif
      i2=i2+i
      if(i2+wztos(nl1).le.slin2)then
      if(slin3.ge.0)then
      j1=stcbs(1)
      call dkar(nl1,j1,jj)
      call cbtx(i2,jj,j1)
      else
      i=stib(xtstrl(0)+26)
      j3=stib(xtstrp(0)+26)
      stxb(i2+1:i2+i)=stcb(j3:j3-1+i)
      i2=i2+i
      endif
      endif
      else
      i=stib(xtstrl(0)+25)
      if(i2+i.le.slin2)then
      j3=stib(xtstrp(0)+25)
      stxb(i2+1:i2+i)=stcb(j3:j3-1+i)
      endif
      i2=i2+i
      if(i2+wztos(nl1).le.slin2)then
      if(slin3.ge.0)then
      j1=stcbs(1)
      call dkar(nl1,j1,jj)
      call cbtx(i2,jj,j1)
      else
      i=stib(xtstrl(0)+26)
      j3=stib(xtstrp(0)+26)
      stxb(i2+1:i2+i)=stcb(j3:j3-1+i)
      i2=i2+i
      endif
      endif
      i=1
      if(i2+i.le.slin2)then
      stxb(i2+1:i2+1)=char(aminus)
      endif
      i2=i2+1
      if(i2+wztos(nl2).le.slin2)then
      if(slin3.ge.0)then
      j1=stcbs(1)
      call dkar(nl2,j1,jj)
      call cbtx(i2,jj,j1)
      else
      i=stib(xtstrl(0)+26)
      j3=stib(xtstrp(0)+26)
      stxb(i2+1:i2+i)=stcb(j3:j3-1+i)
      i2=i2+i
      endif
      endif
      endif
      endif
      if(i2.ge.slin1)then
      j1=slin1
   60 continue
      if(j1.gt.0)then
      if(ichar(stxb(j1:j1)).ne.aspace)then
      j1=j1-1
      goto 60
      endif
      endif
      if(j1.gt.1)then
      print *,stxb(1:j1-1)
      endif
      print *,stcb(1:stab)//stxb(j1:i2)
      else
      print *,stxb(1:i2)
      endif
   70 continue
      call hrul(2)
      print *,stxb(1:ssrec-1)
      print *,''
      jflag(7)=3
      return
      end
      subroutine messag(istop,nl1,nl2,nf1)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( nfiles=5 )
      common/z5in/dunit,munit,ounit,sunit,funit
      common/z12g/jflag(1:11),mflag(1:24)
      character*(srec) auxlin
      common/z22g/auxlin
      integer xunit(nfiles)
      logical lunit
      call mput(istop,nl1,nl2,nf1)
      if(istop.ne.0)then
      xunit(1)=dunit
      xunit(2)=munit
      xunit(3)=funit
      xunit(4)=sunit
      xunit(5)=ounit
      iout=5
      do 10 i=1,nfiles
      if(xunit(i).ne.0)then
      ios=1
      lunit=.false.
      inquire(unit=xunit(i),opened=lunit,iostat=ios)
      if(ios.ne.0)then
      auxlin(1:srec)='system/filesystem error'
      call mput(istop,0,0,0)
      elseif(lunit)then
      ios=1
      if(i.ne.iout)then
      close(unit=xunit(i),status='keep',iostat=ios)
      else
      close(unit=xunit(iout),status='delete',iostat=ios)
      endif
      if(ios.ne.0)then
      auxlin(1:srec)='system/filesystem error'
      call mput(istop,0,0,0)
      endif
      endif
      endif
   10 continue
      stop
      endif
      return
      end
      subroutine hrul(ii)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sxbuff=2040 )
      character*(sxbuff) stxb
      common/z27g/stxb
      common/z51g/aplus,aminus,azero,anine,aeq
      common/z52g/acf0(0:127),acf1(0:127),lfeed,aspace,dquote,squote
      stxb(1:1)=char(aspace)
      j1=ssrec
      if(ii.eq.1)then
      do i1=2,j1
      stxb(i1:i1)=char(aminus)
      enddo
      elseif(ii.eq.2)then
      do i1=3,j1-1
      stxb(i1:i1)=char(aeq)
      enddo
      stxb(2:2)=char(aspace)
      stxb(j1:j1)=char(aspace)
      endif
      return
      end
      subroutine uput(ind)
      implicit integer(a-z)
      save
      parameter ( srec=81, ssrec=62 )
      parameter ( sxbuff=2040 )
      common/z12g/jflag(1:11),mflag(1:24)
      character*(sxbuff) stxb
      common/z27g/stxb
      common/z52g/acf0(0:127),acf1(0:127),lfeed,aspace,dquote,squote
      if(jflag(7).eq.0)then
      print *,''
      endif
      ii=sxbuff-max(80,ssrec)
      if(ii.gt.0)then
      call hrul(2)
      print *,stxb(1:ssrec-1)
      else
      print *,''
      endif
      if(ind.eq.1)then
      print *,'   error: internal buffer size (sibuff) is too small'
      elseif(ind.eq.2)then
      print *,'   error: internal buffer size (scbuff) is too small'
      elseif(ind.eq.3)then
      print *,'   error: largest available integer seems too small'
      elseif(ind.eq.4)then
      print *,'   error: integer too large (in absolute value)'
      elseif(ind.eq.5)then
      print *,'   error: parameter "sxbuff" is not properly set'
      elseif(ind.eq.6)then
      print *,'   error: parameter "sibuff" is not properly set'
      elseif(ind.eq.7)then
      print *,'   error: parameter "scbuff" is not properly set'
      elseif(ind.eq.8)then
      print *,'   error: parameter "mimal" is not properly set'
      elseif(ind.eq.9)then
      print *,'   error: parameter "mcmal" is not properly set'
      elseif(ind.eq.10)then
      print *,'   error: parameter "srec" is not properly set'
      elseif(ind.eq.11)then
      print *,'   error: parameter "ssrec" is not properly set'
      endif
      if(ii.gt.0)then
      print *,stxb(1:ssrec-1)
      endif
      print *,''
      stop
      return
      end
      subroutine wput(sind,j1,j2)
      implicit integer(a-z)
      save
      parameter ( sibuff=524286, scbuff=131064, mimal=2, mcmal=2 )
      parameter ( srec=81, ssrec=62 )
      character*(srec) auxlin
      common/z22g/auxlin
      common/z28g/wera(0:0),werb(0:0),nwer
      common/z30g/stib(1:sibuff)
      character*(scbuff) stcb
      common/z31g/stcb
      common/z35g/stibs(0:mimal),stcbs(0:mcmal),imal(0:0),cmal(0:0)
      aind=abs(sind)
      ii=0
      if(aind.gt.nwer)then
      ii=1
      elseif(j1.lt.0)then
      ii=1
      elseif(j1.eq.0)then
      if(j2.ne.0)then
      ii=1
      endif
      else
      jj=j2-j1
      if(jj.lt.0)then
      ii=1
      elseif(jj.ge.srec-1)then
      ii=1
      endif
      endif
      if(ii.ne.0)then
      auxlin(1:srec)='wput_1'
      call messag(1,0,0,0)
      endif
      i1=stib(wera(0)+aind)
      i2=stib(werb(0)+aind)
      ll=i2-i1
      ii=0
      if(i1.le.0)then
      ii=1
      elseif(i1.ge.stcbs(1))then
      ii=1
      elseif(ll.lt.0)then
      ii=1
      elseif(ll.ge.srec)then
      ii=1
      endif
      if(ii.ne.0)then
      auxlin(1:srec)='wput_2'
      call messag(1,0,0,0)
      endif
      ll=ll+1
      if(j1.ne.0)then
      jj=(ll-j1)+(j2+1)
      else
      jj=0
      endif
      kk=srec-6
      if(jj.gt.0)then
      if(jj.lt.srec)then
      auxlin(1:srec)=stcb(i1:i2)//stcb(j1:j2)
      elseif(ll.le.kk)then
      auxlin(1:srec)=stcb(i1:i2)//stcb(j1:j1+(kk-ll))//' ...'
      else
      auxlin(1:srec)=stcb(i1:i1+kk)//' ...'
      endif
      else
      if(ll.lt.srec)then
      auxlin(1:srec)=stcb(i1:i2)
      else
      auxlin(1:srec)=stcb(i1:i1+kk)//' ...'
      endif
      endif
      if(sind.gt.0)then
      call messag(0,0,0,0)
      else
      call messag(1,0,0,0)
      endif
      return
      end

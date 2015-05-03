C Output from Public domain Ratfor, version 1.0
      subroutine distfun(x,n,p,n2,d)
      implicit double precision (a-h,o-z)
      integer n,p,n2
      real x(n,p),d(n2,p)
      ii=0
      do23000 i=1,n-1
      do23002 ip=i+1,n
      ii=ii+1
      do23004 j=1,p
      d(ii,j)=abs(x(i,j)-x(ip,j))
23004 continue
23005 continue
23002 continue
23003 continue
23000 continue
23001 continue
      return
      end

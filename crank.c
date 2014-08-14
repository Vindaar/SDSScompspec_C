/*************************************************************************** 
CRANK: This is the crank algorithm taken from Numercial Recipes. Here's
what NR has to say about it:

"Given a sorted array w[1..n], replaces the elements by their rank,
including midranking of ties, and returns as s the sum of f3-f, where f is
the number of elements in each tie."

I have changed the raw routine to include the following features:
0) Double precision is used.
1) Function returns an integer.

****************************************************************************/

int crank(unsigned long n, double *w, double *s) {

  unsigned long j=0,ji=0,jt=0;
  double        t=0.0,rank=0.0;

  *s=0.0;
  while (j<n-1) {
    if (w[j+1]!=w[j]) {
      w[j]=(double)(j+1);
      ++j;
    }
    else {
      for (jt=j+1; jt<n && w[jt]==w[j]; jt++);
      rank=0.5*(double)(j+jt+1);
      for (ji=j; ji<=(jt-1); ji++) w[ji]=rank;
      t=(double)(jt-j);
      *s+=t*t*t-t;
      j=jt;
    }
  }
  if (j==n-1) w[n-1]=(double)n;

  return 1;

}

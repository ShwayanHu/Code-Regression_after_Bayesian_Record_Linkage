#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <math.h>
  
void RandomPickNoSort(int *n, double *p, int *ans) {
  double rU;
  int i;
  int nm1 = n[0] - 1;
  
  for (i = 1; i < n[0]; i++)
    p[i] += p[i - 1];
  
  GetRNGstate();
  rU = runif(0, p[nm1]);
  for (i = 0; i < nm1; i++) {
    if (rU <= p[i])
      break;
  }
  ans[0] = i;
  PutRNGstate();
}

void BayesianRecordLinkageWithSeed(
    int *nIter, int *n1, int *n2,
    int *nFields, int *nDisagLevs,
    int *uCompPatts, int *nCompPatts, 
    int *codesCompPatts, int *freqAgrLevs,
    double *a, double *b, 
    double *aBM, double *bBM,
    int *Z_known,  // Input: known links (R-style, 1-based, -1 for unknown)
    int *Z,        // Output: Gibbs samples of links (R-style, 1-based)
    double *m, double *u
) {   
  int iter, i, j, k, f, al, al1, co;
  int from[1];
  int pickedPosition[1];
  int nValidLabels[1];
  int nBinAgrLevels[1];
  int Znew[n2[0]];
  int freeLinks[n1[0]];
  int nFreeLinks[1];
  int freqCodesLinks[nCompPatts[0]];
  double uLambda[nCompPatts[0]];
  double smm[nFields[0]];
  double smu[nFields[0]];
  
  /* Initialize nBinAgrLevels */
  nBinAgrLevels[0] = 0;
  for (f = 0; f < nFields[0]; f++) {
    nBinAgrLevels[0] += nDisagLevs[f];
  }
  
  double logmMINUSlogu[nBinAgrLevels[0]];    
  double A[nBinAgrLevels[0]];
  double B[nBinAgrLevels[0]];
  double mnew[nBinAgrLevels[0]];
  double unew[nBinAgrLevels[0]];
  
  /* Initialize Znew and handle known links */
  for (j = 0; j < n2[0]; j++) {
    if (Z_known[j] != -1) {
      Znew[j] = Z_known[j] - 1;  // Convert to 0-based
    } else {
      Znew[j] = n1[0] + j;       // Unlinked: n1 + j (0-based)
    }
  }
  
  /* Initialize freeLinks (1=free, 0=occupied) */
  for (i = 0; i < n1[0]; i++) {
    freeLinks[i] = 1;
  }
  nFreeLinks[0] = n1[0];
  
  /* Mark known links as occupied */
  for (j = 0; j < n2[0]; j++) {
    if (Z_known[j] != -1) {
      freeLinks[Z_known[j] - 1] = 0;
      nFreeLinks[0]--;
    }
  }
  
  /* Initialize freqCodesLinks */
  for (co = 0; co < nCompPatts[0]; co++) {
    freqCodesLinks[co] = 0;
  }
  
  /* Count known links in freqCodesLinks */
  from[0] = 0;
  for (j = 0; j < n2[0]; j++) {
    if (Z_known[j] != -1) {
      freqCodesLinks[codesCompPatts[from[0] + (Z_known[j] - 1)]]++;
    }
    from[0] += n1[0];
  }
  
  /* Gibbs iterations */
  for (iter = 0; iter < nIter[0]; iter++) {
    /* Compute A and B for m and u */
    for (al = 0; al < nBinAgrLevels[0]; al++) {
      A[al] = 0;
      for (co = 0; co < nCompPatts[0]; co++) {
        if (uCompPatts[nBinAgrLevels[0] * co + al]) 
          A[al] += freqCodesLinks[co];
      }
    }
    
    for (al = 0; al < nBinAgrLevels[0]; al++) {
      B[al] = freqAgrLevs[al] - A[al];
    }
    
    /* Sample mnew and unew from Gamma distributions */
    for (al = 0; al < nBinAgrLevels[0]; al++) {
      mnew[al] = rgamma(A[al] + a[al], 1);
      unew[al] = rgamma(B[al] + b[al], 1);
    }
    
    /* Normalize mnew and unew within each field */
    al1 = 0;
    for (f = 0; f < nFields[0]; f++) {
      smm[f] = 0;
      smu[f] = 0;
      for (al = 0; al < nDisagLevs[f]; al++) {
        smm[f] += mnew[al1];
        smu[f] += unew[al1];
        al1++;
      }
    }
    
    al1 = 0;
    for (f = 0; f < nFields[0]; f++) {
      for (al = 0; al < nDisagLevs[f]; al++) {
        mnew[al1] /= smm[f];
        unew[al1] /= smu[f];
        m[iter * nBinAgrLevels[0] + al1] = mnew[al1];
        u[iter * nBinAgrLevels[0] + al1] = unew[al1];
        al1++;
      }
    }
    
    /* Compute log(m) - log(u) for linkage probabilities */
    for (al = 0; al < nBinAgrLevels[0]; al++) {
      logmMINUSlogu[al] = log(mnew[al]) - log(unew[al]);
    }
    
    /* Compute uLambda for each comparison pattern */
    for (co = 0; co < nCompPatts[0]; co++) {
      uLambda[co] = 0;
      for (al = 0; al < nBinAgrLevels[0]; al++) {
        if (uCompPatts[nBinAgrLevels[0] * co + al]) 
          uLambda[co] += logmMINUSlogu[al];
      }
    }
    
    /* Sample new links (only for unknown records) */
    from[0] = 0;
    for (j = 0; j < n2[0]; j++) {
      /* Skip known links */
      if (Z_known[j] != -1) {
        Z[iter * n2[0] + j] = Z_known[j];  // Known links are already 1-based
        from[0] += n1[0];
        continue;
      }
      
      /* Release old link (if any) */
      if (Znew[j] < n1[0]) {
        freeLinks[Znew[j]] = 1;
        nFreeLinks[0]++;
        freqCodesLinks[codesCompPatts[from[0] + Znew[j]]]--;
      }
      
      /* Sample new link */
      if (nFreeLinks[0] == 0) {
        /* No free links: remain unlinked */
        Znew[j] = n1[0] + j;
      } else {
        /* Compute probabilities for free links */
        double Lambda_i_free[nFreeLinks[0] + 1];
        int which_freelinks[nFreeLinks[0] + 1];
        
        k = 1;
        for (i = 0; i < n1[0]; i++) {
          if (freeLinks[i]) {
            Lambda_i_free[k] = exp(uLambda[codesCompPatts[from[0] + i]]);
            which_freelinks[k] = i;
            k++;
          }
        }
        
        /* Add "no link" option */
        Lambda_i_free[0] = nFreeLinks[0] * (n2[0] - n1[0] + nFreeLinks[0] - 1 + bBM[0]) / 
        (n1[0] - nFreeLinks[0] + aBM[0]);
        nValidLabels[0] = nFreeLinks[0] + 1;
        which_freelinks[0] = n1[0] + j;
        
        /* Randomly pick a link */
        RandomPickNoSort(nValidLabels, Lambda_i_free, pickedPosition);
        Znew[j] = which_freelinks[pickedPosition[0]];
        
        /* Update if linked */
        if (Znew[j] < n1[0]) {
          freeLinks[Znew[j]] = 0;
          nFreeLinks[0]--;
          freqCodesLinks[codesCompPatts[from[0] + Znew[j]]]++;
        }
      }
      
      /* Store Z (convert to R-style: 1-based) */
      Z[iter * n2[0] + j] = Znew[j] + 1;  // n1+j becomes n1+j+1 (e.g., 500 -> 501)
      from[0] += n1[0];
    }
  }
}
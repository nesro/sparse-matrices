sparse-matrices
============================

EN: Impact of sparse matrix storage format on efficiency of multiplication of sparse matrices
CZ: Vliv formátu uložení řídké matice na výkonnost násobení řídkých matic

CZ instructions:
 * Nastudujte formáty uložení řídkých matic COO, CSR, BSR, Quadtree a algoritmy
   pro násobení matic v těchto formátech.
 
 * Navrhněte modifikaci formátu Quadtree snížením výšky stromu, listy budou
   tvořeny podmaticemi ve formátu husté matice nebo ve formátech COO, CSR, BSR.

 * V jazyce C implementujte algoritmy násobení matice vektorem a matice maticí
   ve formátech COO, CSR, BSR, Quadtree, modifikovaný formát Quadtree.

 * Odvoďte jejich časové a paměťové složitosti, porovnejte tyto teoretické
   předpoklady s výkonností a paměťovými nároky jednotlivých implementací.
   
features
============================

 * dense matrix multiplication:
   * naive
   * recursive
   * strassen

 * quadtree matrix multiplication
   * csr leaves optimized
   

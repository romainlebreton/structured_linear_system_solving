<TeXmacs|1.99.5>

<style|generic>

<\body>
  <\hide-preamble>
    <assign|va|<macro|<math|<math-ss|a>>>><assign|vb|<macro|<math|<math-ss|b>>>><assign|vc|<macro|<math|<math-ss|c>>>><assign|ve|<macro|<math|<math-ss|e>>>><assign|vf|<macro|<math|<math-ss|f>>>><assign|vg|<macro|<math|<math-ss|g>>>><assign|vh|<macro|<math|<math-ss|h>>>><assign|vi|<macro|<math|<math-ss|i>>>><assign|vr|<macro|<math|<math-ss|r>>>><assign|vu|<macro|<math|<math-ss|u>>>><assign|vv|<macro|<math|<math-ss|v>>>><assign|vx|<macro|<math|<math-ss|x>>>><assign|vy|<macro|<math|<math-ss|y>>>><assign|mA|<macro|<math|<math-ss|A>>>><assign|mB|<macro|<math|<math-ss|B>>>><assign|mC|<macro|<math|<math-ss|C>>>><assign|mD|<macro|<math|<math-ss|D>>>><assign|mG|<macro|<math|<math-ss|G>>>><assign|mH|<macro|<math|<math-ss|H>>>><assign|mI|<macro|<math|<math-ss|I>>>><assign|mJ|<macro|<math|<math-ss|J>>>><assign|mM|<macro|<math|<math-ss|M>>>><assign|mN|<macro|<math|<math-ss|N>>>><assign|mP|<macro|<math|<math-ss|P>>>><assign|mR|<macro|<math|<math-ss|R>>>><assign|mS|<macro|<math|<math-ss|S>>>><assign|mT|<macro|<math|<math-ss|T>>>><assign|mV|<macro|<math|<math-ss|V>>>><assign|mW|<macro|<math|<math-ss|W>>>><assign|mX|<macro|<math|<math-ss|X>>>><assign|mY|<macro|<math|<math-ss|Y>>>><assign|mZ|<macro|<math|<math-ss|Z>>>><assign|K|<macro|<math|<with|math-font|Bbb|K>>>><assign|Q|<macro|<math|<with|math-font|Bbb|Q>>>><assign|Z|<macro|<math|<with|math-font|Bbb|Z>>>>
  </hide-preamble>

  If the dimension of the Kernel of <math|<mT>> is not one, I propose to add
  a Cauchy matrix below <math|<mA>> instead of adding something to
  <math|<mT>>. We are still adding some constraints to ensure that
  <math|<vx>> is unique and <math|<vx>> relains a solution.

  Assume <math|<mA>=<bmatrix|<tformat|<table|<row|<cell|<mA><rsub|r>>|<cell|<mB>>>>>>>
  with <math|det<around*|(|<mA><rsub|r>|)>\<neq\>0> and consider

  <\equation*>
    <mA><rprime|'>=<bmatrix|<tformat|<cwith|2|2|1|1|cell-bborder|0ln>|<cwith|1|-1|1|1|cell-lborder|0ln>|<cwith|1|-1|1|1|cell-rborder|1ln>|<cwith|1|-1|2|2|cell-lborder|1ln>|<cwith|1|1|1|-1|cell-tborder|0ln>|<cwith|1|1|1|-1|cell-bborder|1ln>|<cwith|2|2|1|-1|cell-tborder|1ln>|<cwith|1|1|1|1|cell-lborder|0ln>|<cwith|1|1|2|2|cell-rborder|0ln>|<table|<row|<cell|<mA><rsub|r>>|<cell|<mB>>>|<row|<cell|<around*|(|<bmatrix|<tformat|<table|<row|<cell|c<rsub|1>>>|<row|<cell|\<vdots\>>>|<row|<cell|c<rsub|s>>>>>>*<bmatrix|<tformat|<table|<row|<cell|l<rsub|1>>|<cell|\<ldots\>>|<cell|l<rsub|r>>>>>>|)>\<odot\><mC><rprime|'>>|<cell|<around*|(|<bmatrix|<tformat|<table|<row|<cell|c<rsub|1>>>|<row|<cell|\<vdots\>>>|<row|<cell|c<rsub|s>>>>>>*<bmatrix|<tformat|<table|<row|<cell|l<rsub|r+1>>|<cell|\<ldots\>>|<cell|l<rsub|r+s>>>>>>|)>\<odot\><mC>>>>>>
  </equation*>

  where <math|c<rsub|1>,\<ldots\>,c<rsub|s>,l<rsub|1>,\<ldots\>,l<rsub|r+s>>
  are indeterminates, and <math|<mC>,<mC><rprime|'>> are Cauchy matrices for
  good vectors <math|<vu>,<vv><rsub|0>> and <math|<vu>,<vv><rsub|1>>
  (<math|<vu>> can be chosen freely, while
  <math|<vv><rsub|0>>,<math|<vv><rsub|1>> depend on the Cauchy structure of
  respectively <math|<mA><rsub|r>> and <math|<mB>>).

  The principal leading minor <math|<mA><rsub|r+1>> is a non zero polynomial
  in <math|k<around*|[|c<rsub|1>,\<ldots\>,c<rsub|s>,l<rsub|1>,\<ldots\>,l<rsub|r+s>|]>>
  because its term in <math|c<rsub|1>*l<rsub|r+1>> is
  <math|c<rsub|1>*l<rsub|r+1>*det<around*|(|<mA><rsub|r>|)>*det<around*|(|<mC><rsub|1>|)>>
  which is non zero since the principal leading minor <math|<mC><rsub|1>> is
  non zero.\ 

  Now the principal leading minor <math|<mA><rsub|r+2>> is a non zero
  polynomial because its term in <math|c<rsub|1>*c<rsub|2>*l<rsub|r+1>*l<rsub|r+2>>
  is <math|c<rsub|1>*c<rsub|2>*l<rsub|r+1>*l<rsub|r+2>*det<around*|(|<mA><rsub|r>|)>*det<around*|(|<mC><rsub|2>|)>\<neq\>0>
  since <math|det<around*|(|<mC><rsub|2>|)>\<neq\>0>. Indeed the presence of
  <math|l<rsub|r+1>*l<rsub|r+2>> in <math|det<around*|(|<mA><rsub|r+2>|)>>
  means that we chose a permutation that involves two terms of the bottom
  right part of <math|<mA><rprime|'>>, and the other terms must come from
  <math|<mA><rsub|r>>.

  \;

  \;
</body>

<initial|<\collection>
</collection>>
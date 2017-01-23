<TeXmacs|1.99.5>

<style|generic>

<\body>
  <\hide-preamble>
    <assign|va|<macro|<math|<math-ss|a>>>>

    <assign|vb|<macro|<math|<math-ss|b>>>>

    <assign|vc|<macro|<math|<math-ss|c>>>>

    <assign|ve|<macro|<math|<math-ss|e>>>>

    <assign|vf|<macro|<math|<math-ss|f>>>>

    <assign|vg|<macro|<math|<math-ss|g>>>>

    <assign|vh|<macro|<math|<math-ss|h>>>>

    <assign|vr|<macro|<math|<math-ss|r>>>>

    <assign|vu|<macro|<math|<math-ss|u>>>>

    <assign|vv|<macro|<math|<math-ss|v>>>>

    <assign|vx|<macro|<math|<math-ss|x>>>>

    <assign|vy|<macro|<math|<math-ss|y>>>>

    <assign|mA|<macro|<math|<math-ss|A>>>>

    <assign|mB|<macro|<math|<math-ss|B>>>>

    <assign|mC|<macro|<math|<math-ss|C>>>>

    <assign|mD|<macro|<math|<math-ss|D>>>>

    <assign|mG|<macro|<math|<math-ss|G>>>>

    <assign|mH|<macro|<math|<math-ss|H>>>>

    <assign|mJ|<macro|<math|<math-ss|J>>>>

    <assign|mM|<macro|<math|<math-ss|M>>>>

    <assign|mN|<macro|<math|<math-ss|N>>>>

    <assign|mS|<macro|<math|<math-ss|S>>>>

    <assign|mT|<macro|<math|<math-ss|T>>>>

    <assign|mV|<macro|<math|<math-ss|V>>>>

    <assign|mW|<macro|<math|<math-ss|W>>>>

    <assign|mX|<macro|<math|<math-ss|X>>>>

    <assign|mY|<macro|<math|<math-ss|Y>>>>

    <assign|mZ|<macro|<math|<math-ss|Z>>>>

    <assign|K|<macro|<math|<with|math-font|Bbb|K>>>>

    <assign|Q|<macro|<math|<with|math-font|Bbb|Q>>>>

    <assign|Z|<macro|<math|<with|math-font|Bbb|Z>>>>

    <assign|M|<macro|<math|<with|math-font|cal*|M>>>>
  </hide-preamble>

  <assign|mA|<macro|<math|<math-ss|A>>>><section*|From Toeplitz to Cauchy
  matrices>

  Let <math|<mZ>=<matrix|<tformat|<table|<row|<cell|0>|<cell|0>|<cell|0>|<cell|0>>|<row|<cell|1>|<cell|\<ddots\>>|<cell|\<ddots\>>|<cell|0>>|<row|<cell|0>|<cell|\<ddots\>>|<cell|\<ddots\>>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|1>|<cell|0>>>>>>
  be the displacement matrix that acts as

  <\equation*>
    <mZ>*<mA>=<mA>\<downarrow\>,<htab|5mm><mZ>*<rsup|t><mA>=<mA>\<uparrow\>,<htab|5mm><mA>*<mZ>=<mA>\<leftarrow\>,<htab|5mm><mA>*<mZ><rsup|t>=<mA>\<rightarrow\>.
  </equation*>

  Toeplitz and Cauchy matrices are introduced in the paper. To fix the
  notations, denote <math|<mV><rsub|<vu>>>,
  <math|<wide|<mV>|\<bar\>><rsub|<vu>>>, <math|<mW><rsub|<vu>>> and
  <math|<wide|<mW>|\<bar\>><rsub|<vu>>> be the following \PVandermonde\Q
  matrices

  <\equation*>
    <mV><rsub|<vu>>=<matrix|<tformat|<table|<row|<cell|1>|<cell|u<rsub|1>>|<cell|\<cdots\>>|<cell|u<rsub|1><rsup|n-1>>>|<row|<cell|\<vdots\>>|<cell|\<vdots\>>|<cell|>|<cell|\<vdots\>>>|<row|<cell|1>|<cell|u<rsub|n>>|<cell|\<cdots\>>|<cell|u<rsub|n><rsup|n-1>>>>>>,<htab|5mm><wide|<mV>|\<bar\>><rsub|<vu>>=<matrix|<tformat|<table|<row|<cell|u<rsub|1><rsup|n-1>>|<cell|*\<cdots\>>|<cell|u<rsub|1>>|<cell|1>>|<row|<cell|\<vdots\>>|<cell|>|<cell|\<vdots\>>|<cell|\<vdots\>*>>|<row|<cell|u<rsub|n><rsup|n-1>>|<cell|\<cdots\>>|<cell|u<rsub|n>>|<cell|1>>>>>,
  </equation*>

  <\equation*>
    <mW><rsub|<vu>>=<matrix|<tformat|<table|<row|<cell|u<rsub|1><rsup|n-1>>|<cell|*\<cdots\>>|<cell|u<rsub|n><rsup|n-1>>>|<row|<cell|\<vdots\>>|<cell|>|<cell|\<vdots\>>>|<row|<cell|u<rsub|1>>|<cell|\<cdots\>>|<cell|u<rsub|n>>>|<row|<cell|1>|<cell|\<cdots\>>|<cell|1>>>>>,<htab|5mm><wide|<mW>|\<bar\>><rsub|<vu>>=<matrix|<tformat|<table|<row|<cell|1>|<cell|\<cdots\>>|<cell|1>>|<row|<cell|u<rsub|1>>|<cell|\<cdots\>>|<cell|u<rsub|n>>>|<row|<cell|\<vdots\>>|<cell|>|<cell|\<vdots\>>>|<row|<cell|u<rsub|1><rsup|n-1>>|<cell|*\<cdots\>>|<cell|u<rsub|n><rsup|n-1>>>>>>.
  </equation*>

  I keep the paper notations for <math|<mV><rsub|<vu>>>,
  <math|<mW><rsub|<vu>>> and introduce <math|<wide|<mV>|\<bar\>><rsub|<vu>>>,
  <math|<wide|<mW>|\<bar\>><rsub|<vu>>> which are the correct matrices in my
  opinion (see below).

  We have the following relations

  <\eqnarray*>
    <tformat|<table|<row|<cell|<mD><rsub|<vu>>*<mV><rsub|<vu>>>|<cell|=>|<cell|<mV><rsub|<vu>>*<mZ>+<matrix|<tformat|<table|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>|<cell|u<rsub|1><rsup|n>>>|<row|<cell|\<vdots\>>|<cell|>|<cell|\<vdots\>>|<cell|\<vdots\>>>|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>|<cell|u<rsub|n><rsup|n>>>>>>>>|<row|<cell|<mD><rsub|<vu>>*<wide|<mV>|\<bar\>><rsub|<vu>>>|<cell|=>|<cell|<wide|<mV>|\<bar\>><rsub|<vu>>*<mZ><rsup|t>+<matrix|<tformat|<table|<row|<cell|u<rsub|1><rsup|n>>|<cell|0>|<cell|\<cdots\>>|<cell|0>>|<row|<cell|\<vdots\>>|<cell|\<vdots\>>|<cell|>|<cell|\<vdots\>>>|<row|<cell|u<rsub|n><rsup|n>>|<cell|0>|<cell|\<cdots\>>|<cell|0>>>>>>>|<row|<cell|<mW><rsub|<vu>>*<mD><rsub|<vu>>>|<cell|=>|<cell|<mZ>*<mW><rsub|<vu>>+<matrix|<tformat|<table|<row|<cell|u<rsub|1><rsup|n>>|<cell|*\<cdots\>>|<cell|u<rsub|n><rsup|n>>>|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>>|<row|<cell|\<vdots\>>|<cell|>|<cell|\<vdots\>>>|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>>>>>>>|<row|<cell|<wide|<mW>|\<bar\>><rsub|<vu>>*<mD><rsub|<vu>>>|<cell|=>|<cell|<mZ><rsup|t>*<wide|<mW>|\<bar\>><rsub|<vu>>+<matrix|<tformat|<table|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>>|<row|<cell|\<vdots\>>|<cell|>|<cell|\<vdots\>>>|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>>|<row|<cell|u<rsub|1><rsup|n>>|<cell|*\<cdots\>>|<cell|u<rsub|n><rsup|n>>>>>>>>>>
  </eqnarray*>

  \;

  If <math|<mA><rprime|'>=<wide|<mV>|\<bar\>><rsub|<vu>>*<mA>*<wide|<mW>|\<bar\>><rsub|<vv>>>
  then

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<nabla\><rsub|<vu>,<vv>><around*|(|<mA><rprime|'>|)>>|<cell|=>|<cell|<mD><rsub|<vu>>*<wide|<mV>|\<bar\>><rsub|<vu>>*<mA>*<wide|<mW>|\<bar\>><rsub|<vv>>-<wide|<mV>|\<bar\>><rsub|<vu>>*<mA>*<wide|<mW>|\<bar\>><rsub|<vv>>*<mD><rsub|<vv>>>>|<row|<cell|>|<cell|=>|<cell|<wide|<mV>|\<bar\>><rsub|<vu>>*<around*|(|<mZ><rsup|t>*<mA>-<mA>*<mZ><rsup|t>|)>*<wide|<mW>|\<bar\>><rsub|<vv>>+<matrix|<tformat|<table|<row|<cell|u<rsub|1><rsup|n>>|<cell|0>|<cell|\<cdots\>>|<cell|0>>|<row|<cell|\<vdots\>>|<cell|\<vdots\>>|<cell|>|<cell|\<vdots\>>>|<row|<cell|u<rsub|n><rsup|n>>|<cell|0>|<cell|\<cdots\>>|<cell|0>>>>>*<mA>*<wide|<mW>|\<bar\>><rsub|<vv>>-<wide|<mV>|\<bar\>><rsub|<vu>>*<mA>*<matrix|<tformat|<table|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>>|<row|<cell|\<vdots\>>|<cell|>|<cell|\<vdots\>>>|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>>|<row|<cell|v<rsub|1><rsup|n>>|<cell|*\<cdots\>>|<cell|v<rsub|n><rsup|n>>>>>>.>>>>
  </eqnarray*>

  Otherwise ff <math|<mA><rprime|'>=<mV><rsub|<vu>>*<mA>*<mW><rsub|<vv>>>
  then

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<nabla\><rsub|<vu>,<vv>><around*|(|<mA><rprime|'>|)>>|<cell|=>|<cell|<mD><rsub|<vu>>*<mV><rsub|<vu>>*<mA>*<mW><rsub|<vv>>-<mV><rsub|<vu>>*<mA>*<mW><rsub|<vv>>*<mD><rsub|<vv>>>>|<row|<cell|>|<cell|=>|<cell|<mV><rsub|<vu>>*<around*|(|<mZ>*<mA>-<mA>*<mZ>|)>*<mW><rsub|<vv>>+<matrix|<tformat|<table|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>|<cell|u<rsub|1><rsup|n>>>|<row|<cell|\<vdots\>>|<cell|>|<cell|\<vdots\>>|<cell|\<vdots\>>>|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>|<cell|u<rsub|n><rsup|n>>>>>>*<mA>*<mW><rsub|<vv>>-<mV><rsub|<vu>>*<mA><matrix|<tformat|<table|<row|<cell|u<rsub|1><rsup|n>>|<cell|*\<cdots\>>|<cell|u<rsub|n><rsup|n>>>|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>>|<row|<cell|\<vdots\>>|<cell|>|<cell|\<vdots\>>>|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>>>>>.>>>>
  </eqnarray*>

  Or we can express it using the more common
  <math|<around*|(|<mA>-<mZ>*<mA>*<mZ><rsup|t>|)>> displacement operator
  using

  <\equation*>
    <wide|<mV>|\<bar\>><rsub|<vu>>=<mD><rsub|<vu>>*<wide|<mV>|\<bar\>><rsub|<vu>>*<mZ>+<matrix|<tformat|<table|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>|<cell|1>>|<row|<cell|\<vdots\>>|<cell|>|<cell|\<vdots\>>|<cell|\<vdots\>>>|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>|<cell|1>>>>>
  </equation*>

  to get\ 

  <\eqnarray*>
    <tformat|<table|<row|<cell|\<nabla\><rsub|<vu>,<vv>><around*|(|<mA><rprime|'>|)>>|<cell|=>|<cell|<mD><rsub|<vu>>*<wide|<mV>|\<bar\>><rsub|<vu>>*<mA>*<wide|<mW>|\<bar\>><rsub|<vv>>-<wide|<mV>|\<bar\>><rsub|<vu>>*<mA>*<wide|<mW>|\<bar\>><rsub|<vv>>*<mD><rsub|<vv>>>>|<row|<cell|>|<cell|=>|<cell|<mD><rsub|<vu>>*<wide|<mV>|\<bar\>><rsub|<vu>>*<mA>*<wide|<mW>|\<bar\>><rsub|<vv>>-<around*|[|<mD><rsub|<vu>>*<wide|<mV>|\<bar\>><rsub|<vu>>*<mZ>+<matrix|<tformat|<table|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>|<cell|1>>|<row|<cell|\<vdots\>>|<cell|>|<cell|\<vdots\>>|<cell|\<vdots\>>>|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>|<cell|1>>>>>|]>*<mA>*<wide|<mW>|\<bar\>><rsub|<vv>>*<mD><rsub|<vv>>>>|<row|<cell|>|<cell|=>|<cell|<mD><rsub|<vu>>*<wide|<mV>|\<bar\>><rsub|<vu>>*<mA>*<wide|<mW>|\<bar\>><rsub|<vv>>-<mD><rsub|<vu>>*<wide|<mV>|\<bar\>><rsub|<vu>>*<mZ>*<mA>*<wide|<mW>|\<bar\>><rsub|<vv>>*<mD><rsub|<vv>>-<matrix|<tformat|<table|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>|<cell|1>>|<row|<cell|\<vdots\>>|<cell|>|<cell|\<vdots\>>|<cell|\<vdots\>>>|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>|<cell|1>>>>>*<mA>*<wide|<mW>|\<bar\>><rsub|<vv>>*<mD><rsub|<vv>>>>|<row|<cell|>|<cell|=>|<cell|<mD><rsub|<vu>>*<wide|<mV>|\<bar\>><rsub|<vu>>*<mA>*<wide|<mW>|\<bar\>><rsub|<vv>>-<mD><rsub|<vu>>*<wide|<mV>|\<bar\>><rsub|<vu>>*<mZ>*<mA>*<around*|[|<mZ><rsup|t>*<wide|<mW>|\<bar\>><rsub|<vv>>+<matrix|<tformat|<table|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>>|<row|<cell|\<vdots\>>|<cell|>|<cell|\<vdots\>>>|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>>|<row|<cell|v<rsub|1><rsup|n>>|<cell|*\<cdots\>>|<cell|v<rsub|n><rsup|n>>>>>>|]>-<matrix|<tformat|<table|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>|<cell|1>>|<row|<cell|\<vdots\>>|<cell|>|<cell|\<vdots\>>|<cell|\<vdots\>>>|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>|<cell|1>>>>>*<mA>*<wide|<mW>|\<bar\>><rsub|<vv>>*<mD><rsub|<vv>>>>|<row|<cell|>|<cell|=>|<cell|<mD><rsub|<vu>>*<wide|<mV>|\<bar\>><rsub|<vu>>*<around*|(|<mA>-<mZ>*<mA>*<mZ><rsup|t>|)>*<wide|<mW>|\<bar\>><rsub|<vv>>-<mD><rsub|<vu>>*<wide|<mV>|\<bar\>><rsub|<vu>>*<mZ>*<mA>*<matrix|<tformat|<table|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>>|<row|<cell|\<vdots\>>|<cell|>|<cell|\<vdots\>>>|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>>|<row|<cell|v<rsub|1><rsup|n>>|<cell|*\<cdots\>>|<cell|v<rsub|n><rsup|n>>>>>>-<matrix|<tformat|<table|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>|<cell|1>>|<row|<cell|\<vdots\>>|<cell|>|<cell|\<vdots\>>|<cell|\<vdots\>>>|<row|<cell|0>|<cell|\<cdots\>>|<cell|0>|<cell|1>>>>>*<mA>*<wide|<mW>|\<bar\>><rsub|<vv>>*<mD><rsub|<vv>>.>>>>
  </eqnarray*>

  \;

  <chapter*|Old stuff>

  Displacement operators :

  <\enumerate>
    <item>Toeplitz:\ 

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<phi\><rsup|+>>|<cell|:>|<cell|<mA>\<mapsto\><mA>-<mZ>*<mA>*<mZ><rsup|t>=<mA>-<around|(|<mA><text|<nbsp>shifted
      down and right by one place>|)>>>|<row|<cell|\<phi\>>|<cell|:>|<cell|<mA>\<mapsto\><mZ><rsup|t>*<mA>-<mA>*<mZ><rsup|t>=<around*|(|<mA><text|<nbsp>shifted
      up by one place>|)>-<around|(|<mA><text|<nbsp>shifted right by one
      place>|)>>>>>
    </eqnarray*>

    <item>Cauchy:

    <\equation*>
      \<nabla\><rsub|<vu>,<vv>>:<mA>\<mapsto\><mD><rsub|<vu>>*<mA>-<mA>*<mD><rsub|<vv>>
    </equation*>

    <item>Vandermonde:

    <\eqnarray*>
      <tformat|<table|<row|<cell|\<varphi\><rsup|+><rsub|<vu>>>|<cell|:>|<cell|<mA>\<mapsto\><mD><rsub|<vu>>*<mA>-<mA>*<mZ><rsup|t>>>|<row|<cell|\<varphi\><rsub|<vu>>>|<cell|:>|<cell|<mA>\<mapsto\><mZ><rsup|t>*<mA>-<mA>*<mD><rsub|<vu>>>>>>
    </eqnarray*>
  </enumerate>
</body>

<initial|<\collection>
</collection>>

<\references>
  <\collection>
    <associate|auto-1|<tuple|?|?>>
    <associate|auto-2|<tuple|?|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|From
      Toeplitz to Cauchy matrices> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|2fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|font-size|<quote|1.19>|Old
      stuff> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|1fn>
    </associate>
  </collection>
</auxiliary>
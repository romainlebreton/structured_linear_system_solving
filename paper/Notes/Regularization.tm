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

  Let <math|I=<around*|{|1,\<ldots\>,i|}>> for
  <math|i\<leqslant\>rank<around*|(|<mA>|)>>. Then

  <\equation*>
    det<around*|(|<mD><rsub|<vx>>*<mV><rsub|<vu>>*<mA>*<mW><rsub|<vv>>*<mD><rsub|<vy>>|)><rsub|I,I>=<big|sum><rsub|<stack|<tformat|<table|<row|<cell|<around*|\||J|\|>=i>>|<row|<cell|<around*|\||K|\|>=i>>>>>>det<around*|(|<mD><rsub|<vx>>*<mV><rsub|<vu>>|)><rsub|I,J>*det<around*|(|<mA>|)><rsub|J,K>*det<around*|(|<mW><rsub|<vv>>*<mD><rsub|<vy>>|)><rsub|K,I>
  </equation*>

  but <math|det<around*|(|<mD><rsub|<vx>>*<mV><rsub|<vu>>|)><rsub|I,J>=x<rsub|1>*\<cdots\>*x<rsub|i>*det<around*|(|<mV><rsub|<vu>>|)><rsub|I,J>>.
  So the monomial <math|x<rsub|1>*\<cdots\>*x<rsub|i>> do not depend on
  <math|J> and is not unique in the sum.

  In fact, Pan uses the following regularization in his book
  <math|<mA><rprime|''>=<mC><rsub|<va>,<vu>>*<mD><rsub|<vx>>*<mA><rprime|'>*<mD><rsub|<vy>>*<mC><rsub|<vv>,<vb>>>
  since

  <\equation*>
    det<around*|(|<mC><rsub|<va>,<vu>>*<mD><rsub|<vx>>|)><rsub|I,J>=x<rsub|j<rsub|1>>*\<cdots\>*x<rsub|j<rsub|i>>*det<around*|(|<mC><rsub|<va>,<vu>>|)><rsub|I,J>.
  </equation*>
</body>

<initial|<\collection>
</collection>>
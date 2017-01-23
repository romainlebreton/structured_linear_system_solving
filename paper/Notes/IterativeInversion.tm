<TeXmacs|1.99.5>

<style|generic>

<\body>
  <strong|Fact to explain: >It is the same transformations to derive
  <math|S<rsup|<around*|(|i+j|)>>> from <math|S<rsup|<around*|(|j|)>>> than
  to derive <math|S<rsup|<around*|(|i|)>>> from
  <math|S<rsup|<around*|(|0|)>>=A\<in\>\<bbb-K\><rsup|m\<times\>n>> ?

  Let <math|p=min<around*|(|m,n|)>> and consider the matrix

  <\equation*>
    R<rsup|<around*|(|0|)>>=<matrix|<tformat|<table|<row|<cell|A>>|<row|<cell|I<rsub|p,n>>>>>>=<matrix|<tformat|<table|<row|<cell|S<rsup|<around*|(|0|)>>>>|<row|<cell|I<rsub|p,n>>>>>>\<in\>\<bbb-K\><rsup|<around*|(|m+p|)>\<times\>n>
  </equation*>

  where <math|I<rsub|p,n>\<in\>\<bbb-K\><rsup|p\<times\>n>> is the
  rectangular matrix with ones on the diagonal. For
  <math|i\<in\><around*|{|0,\<ldots\>,p|}>>, write <math|A> as

  <\equation*>
    A=<matrix|<tformat|<table|<row|<cell|A<rsub|0,0><rsup|<around*|(|i|)>>>|<cell|A<rsub|0,1><rsup|<around*|(|i|)>>>>|<row|<cell|A<rsub|1,0><rsup|<around*|(|i|)>>>|<cell|A<rsub|1,1><rsup|<around*|(|i|)>>>>>>>
  </equation*>

  where <math|A<rsub|0,0><rsup|<around*|(|i|)>>\<in\>\<bbb-K\><rsup|i\<times\>i>>.\ 

  Now assume that the principal minor <math|A<rsub|0,0><rsup|<around*|(|i|)>>>
  is invertible. Using column operations, we want to transform
  <math|<matrix|<tformat|<table|<row|<cell|A<rsub|0,0><rsup|<around*|(|i|)>>>|<cell|A<rsub|0,1><rsup|<around*|(|i|)>>>>>>>>
  to <math|<matrix|<tformat|<table|<row|<cell|Id<rsub|i>>|<cell|0>>>>>>. When
  we apply this transformation to <math|R<rsup|<around*|(|0|)>>>, we get

  <\equation*>
    R<rsup|<around*|(|i|)>>\<assign\><matrix|<tformat|<table|<row|<cell|Id<rsub|i>>|<cell|0>>|<row|<cell|A<rsub|1,0><rsup|<around*|(|i|)>>*A<rsub|0,0><rsup|<around*|(|i|)>><phantom|><rsup|-1>>|<cell|A<rsub|1,1>-A<rsub|1,0><rsup|<around*|(|i|)>>*A<rsub|0,0><rsup|<around*|(|i|)>><phantom|><rsup|-1>*A<rsub|0,1><rsup|<around*|(|i|)>>>>|<row|<cell|A<rsub|0,0><rsup|<around*|(|i|)>><phantom|><rsup|-1>>|<cell|-A<rsub|0,0><rsup|<around*|(|i|)>><phantom|><rsup|-1>*A<rsub|0,1><rsup|<around*|(|i|)>>>>|<row|<cell|0>|<cell|I<rsub|p-i,n-i>>>>>>=<matrix|<tformat|<table|<row|<cell|A<rsub|0,0><rsup|<around*|(|i|)>>>|<cell|A<rsub|0,1><rsup|<around*|(|i|)>>>>|<row|<cell|A<rsub|1,0><rsup|<around*|(|i|)>>>|<cell|A<rsub|1,1><rsup|<around*|(|i|)>>>>|<row|<cell|Id<rsub|i>>|<cell|0>>|<row|<cell|0>|<cell|I<rsub|p-i,n-i>>>>>>\<times\><matrix|<tformat|<table|<row|<cell|A<rsub|0,0><rsup|<around*|(|i|)>><phantom|><rsup|-1>>|<cell|-A<rsub|0,0><rsup|<around*|(|i|)>><phantom|><rsup|-1>*A<rsub|0,1><rsup|<around*|(|i|)>>>>|<row|<cell|0>|<cell|Id<rsub|n-i>>>>>>
  </equation*>

  Define <math|S<rsup|<around*|(|i|)>>> as the central part of
  <math|R<rsup|<around*|(|i|)>>>, so that

  <\equation*>
    R<rsup|<around*|(|i|)>>=<matrix|<tformat|<cwith|1|1|1|1|cell-tborder|0ln>|<cwith|1|1|1|1|cell-bborder|0ln>|<cwith|2|2|1|1|cell-tborder|0ln>|<cwith|1|1|1|1|cell-lborder|0ln>|<cwith|1|1|1|1|cell-rborder|1ln>|<cwith|1|1|2|2|cell-lborder|1ln>|<cwith|4|4|2|2|cell-tborder|0ln>|<cwith|2|2|2|2|cell-bborder|0ln>|<cwith|4|4|2|2|cell-bborder|0ln>|<cwith|4|4|2|2|cell-lborder|1ln>|<cwith|4|4|1|1|cell-rborder|1ln>|<cwith|4|4|2|2|cell-rborder|0ln>|<cwith|2|2|1|1|cell-row-span|2>|<cwith|2|2|1|1|cell-col-span|2>|<cwith|2|2|1|1|cell-vpart|2>|<cwith|2|2|1|1|cell-valign|c>|<cwith|2|2|1|1|cell-height|2em>|<cwith|2|2|1|1|cell-vmode|exact>|<table|<row|<cell|Id<rsub|i>>|<cell|0>>|<row|<cell|<large|S<rsup|<around*|(|i|)>>>>|<cell|>>|<row|<cell|>|<cell|>>|<row|<cell|0>|<cell|I<rsub|p-i,n-i>>>>>>.
  </equation*>

  <strong|Composition rule.> We want to prove that it is the same
  transformations to derive <math|S<rsup|<around*|(|i+j|)>>> from
  <math|S<rsup|<around*|(|j|)>>> than to derive
  <math|S<rsup|<around*|(|i|)>>> from <math|S<rsup|<around*|(|0|)>>=A\<in\>\<bbb-K\><rsup|m\<times\>n>>.
  Put it differently, we can compute <math|S<rsup|<around*|(|i+j|)>>> as the
  composition of one transformation of order <math|i> and one of order
  <math|j>.

  <\proof>
    The transformation that we apply to <math|R<rsup|<around*|(|0|)>>> to get
    <math|R<rsup|<around*|(|i|)>>> is uniquely determined by the first
    <math|i> rows \ <math|<matrix|<tformat|<table|<row|<cell|A<rsub|0,0><rsup|<around*|(|i|)>>>|<cell|A<rsub|0,1><rsup|<around*|(|i|)>>>>>>>>
    of <math|A>. It consists of a column reduction that maps
    <math|<matrix|<tformat|<table|<row|<cell|A<rsub|0,0><rsup|<around*|(|i|)>>>|<cell|A<rsub|0,1><rsup|<around*|(|i|)>>>>>>>>
    to <math|<matrix|<tformat|<table|<row|<cell|Id<rsub|i>>|<cell|0>>>>>>
    using <math|A<rsub|0,0><rsup|<around*|(|i|)>>> as pivot in the following
    manner

    <\equation*>
      <matrix|<tformat|<table|<row|<cell|A<rsub|0,0><rsup|<around*|(|i|)>>>|<cell|A<rsub|0,1><rsup|<around*|(|i|)>>>>>>>\<times\><matrix|<tformat|<table|<row|<cell|U>|<cell|V>>|<row|<cell|0>|<cell|Id<rsub|n-i>>>>>>=<matrix|<tformat|<table|<row|<cell|Id<rsub|i>>|<cell|0>>>>>,
    </equation*>

    for <math|U=A<rsub|0,0><rsup|<around*|(|i|)>><phantom|><rsup|-1>\<in\>\<bbb-K\><rsup|i\<times\>i>>
    and <math|V=-A<rsub|0,0><rsup|<around*|(|i|)>><phantom|><rsup|-1>*A<rsub|0,1><rsup|<around*|(|i|)>>\<in\>\<bbb-K\><rsup|i\<times\>n-i>>.
    The matrices <math|U> and <math|V> are uniquely determined by last
    equation ; it can be seen easily on the following implied equation

    <\equation*>
      <matrix|<tformat|<table|<row|<cell|A<rsub|0,0><rsup|<around*|(|i|)>>>|<cell|A<rsub|0,1><rsup|<around*|(|i|)>>>>|<row|<cell|0>|<cell|Id<rsub|n-i>>>>>>\<times\><matrix|<tformat|<table|<row|<cell|U>|<cell|V>>|<row|<cell|0>|<cell|Id<rsub|n-i>>>>>>=<matrix|<tformat|<table|<row|<cell|Id<rsub|i>>|<cell|0>>|<row|<cell|0>|<cell|Id<rsub|n-i>>>>>><htab|5mm>\<Rightarrow\><htab|5mm><matrix|<tformat|<table|<row|<cell|U>|<cell|V>>|<row|<cell|0>|<cell|Id<rsub|n-i>>>>>>=<matrix|<tformat|<table|<row|<cell|A<rsub|0,0><rsup|<around*|(|i|)>>>|<cell|A<rsub|0,1><rsup|<around*|(|i|)>>>>|<row|<cell|0>|<cell|Id<rsub|n-i>>>>>><rsup|-1>.
    </equation*>

    To see that this transformation of order <math|i> coincides with the
    composition of column reduction of order <math|k> and
    <math|<around*|(|i-k|)>> for any <math|k\<in\><around*|{|0,\<ldots\>,i|}>>,
    it remains to prove that the composition performs\ 

    <\equation*>
      <matrix|<tformat|<table|<row|<cell|A<rsub|0,0><rsup|<around*|(|i|)>>>|<cell|A<rsub|0,1><rsup|<around*|(|i|)>>>>>>>\<times\><matrix|<tformat|<table|<row|<cell|U>|<cell|V>>|<row|<cell|0>|<cell|Id<rsub|n-i>>>>>>=<matrix|<tformat|<table|<row|<cell|Id<rsub|i>>|<cell|0>>>>>,
    </equation*>

    for some <math|U\<in\>\<bbb-K\><rsup|i\<times\>i>> and
    <math|V\<in\>\<bbb-K\><rsup|i\<times\>n-i>>. Let us write

    <\equation*>
      <matrix|<tformat|<table|<row|<cell|A<rsub|0,0><rsup|<around*|(|i|)>>>|<cell|A<rsub|0,1><rsup|<around*|(|i|)>>>>>>>=<matrix|<tformat|<table|<row|<cell|B<rsub|0,0>>|<cell|B<rsub|0,1>>|<cell|B<rsub|0,2>>>|<row|<cell|B<rsub|1,0>>|<cell|B<rsub|1,1>>|<cell|B<rsub|1,2>>>>>>
    </equation*>

    where <math|B<rsub|0,0>\<in\>\<bbb-K\><rsup|k\<times\>k>> and
    <math|B<rsub|1,1>\<in\>\<bbb-K\><rsup|<around*|(|i-k|)>\<times\><around*|(|i-k|)>>>.
    The first column reduction of order <math|k> performs

    <\equation*>
      <matrix|<tformat|<table|<row|<cell|B<rsub|0,0>>|<cell|B<rsub|0,1>>|<cell|B<rsub|0,2>>>|<row|<cell|B<rsub|1,0>>|<cell|B<rsub|1,1>>|<cell|B<rsub|1,2>>>>>>\<times\><matrix|<tformat|<table|<row|<cell|\<ast\>>|<cell|\<ast\>>|<cell|\<ast\>>>|<row|<cell|0>|<cell|Id<rsub|i-k>>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|Id<rsub|n-i>>>>>>=<matrix|<tformat|<table|<row|<cell|Id<rsub|k>>|<cell|0>|<cell|0>>|<row|<cell|C<rsub|1,0>>|<cell|C<rsub|1,1>>|<cell|C<rsub|1,2>>>>>>.
    </equation*>

    It is followed by the following column reduction of order
    <math|<around*|(|i-k|)>> on the last <math|<around*|(|i-k|)>> rows using
    <math|C<rsub|1,1>> as pivot

    <\equation*>
      <matrix|<tformat|<table|<row|<cell|Id<rsub|k>>|<cell|0>|<cell|0>>|<row|<cell|C<rsub|1,0>>|<cell|C<rsub|1,1>>|<cell|C<rsub|1,2>>>>>>\<times\><matrix|<tformat|<table|<row|<cell|Id<rsub|k>>|<cell|0>|<cell|0>>|<row|<cell|\<ast\>>|<cell|\<ast\>>|<cell|\<ast\>>>|<row|<cell|0>|<cell|0>|<cell|Id<rsub|n-i>>>>>>=<matrix|<tformat|<table|<row|<cell|Id<rsub|k>>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|Id<rsub|i-k>>|<cell|0>>>>>.
    </equation*>

    Since\ 

    <\equation*>
      <matrix|<tformat|<table|<row|<cell|\<ast\>>|<cell|\<ast\>>|<cell|\<ast\>>>|<row|<cell|0>|<cell|Id<rsub|i-k>>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|Id<rsub|n-i>>>>>>\<times\><matrix|<tformat|<table|<row|<cell|Id<rsub|k>>|<cell|0>|<cell|0>>|<row|<cell|\<ast\>>|<cell|\<ast\>>|<cell|\<ast\>>>|<row|<cell|0>|<cell|0>|<cell|Id<rsub|n-i>>>>>>=<matrix|<tformat|<table|<row|<cell|\<ast\>>|<cell|\<ast\>>|<cell|\<ast\>>>|<row|<cell|\<ast\>>|<cell|\<ast\>>|<cell|\<ast\>>>|<row|<cell|0>|<cell|0>|<cell|Id<rsub|n-i>>>>>>
    </equation*>

    the composed operation verifies the same assumption and coincide with the
    transformation of order <math|i> by uniqueness.
  </proof>

  <strong|Link with the inverse. >An important fact with the matrices
  <math|R<rsup|<around*|(|i|)>>> and <math|S<rsup|<around*|(|i|)>>> is that,
  if <math|A\<in\>\<bbb-K\><rsup|n\<times\>n>> is invertible, then

  <\equation*>
    R<rsup|<around*|(|n|)>>=<matrix|<tformat|<table|<row|<cell|Id<rsub|n>>>|<row|<cell|A<rsup|-1>>>>>><math-up|
    and >S<rsup|<around*|(|n|)>>=A<rsup|-1>.
  </equation*>

  More generally, if <math|r=rank<around*|(|A|)>> and
  <math|A<rsub|0,0><rsup|<around*|(|r|)>>> is invertible, then\ 

  <\equation*>
    R<rsup|<around*|(|r|)>>\<assign\><matrix|<tformat|<table|<row|<cell|Id<rsub|i>>|<cell|0>>|<row|<cell|A<rsub|1,0><rsup|<around*|(|r|)>>*A<rsub|0,0><rsup|<around*|(|r|)>><phantom|><rsup|-1>>|<cell|0>>|<row|<cell|A<rsub|0,0><rsup|<around*|(|r|)>><phantom|><rsup|-1>>|<cell|-A<rsub|0,0><rsup|<around*|(|r|)>><phantom|><rsup|-1>*A<rsub|0,1><rsup|<around*|(|r|)>>>>|<row|<cell|0>|<cell|I<rsub|p-r,n-r>>>>>>
  </equation*>

  <strong|Iterative inversion. >The link between the matrices
  <math|R<rsup|<around*|(|i|)>>> and the inverse, together with the
  composition rule, are the keys of Cardinal's paper for Cauchy structured
  matrix inversion. It must also be linked to Strassen's matrix inversion
  formula.
</body>

<initial|<\collection>
</collection>>
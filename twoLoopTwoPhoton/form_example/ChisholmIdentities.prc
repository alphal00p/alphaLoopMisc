#procedure ChisholmIdentities()
  label retry;
    repeat id g(n1?,n2?)*g(n2?,n3?) = g(n1,n3);
    repeat id g(n1?,n2?)*Gtrace(?a,n1?,?b) = Gtrace(?a,n2,?b);
    id once ifmatch->retry Gtrace(?a,n?,n?,?b)=Gtrace(?a,?b)*rat(4-2*ep,1);
    id once ifmatch->retry Gtrace(?a,n?,n1?,n?,?b)=rat(-2+2*ep,1)*Gtrace(?a,n1,?b);
    id once ifmatch->retry Gtrace(?a,n?,n1?,n2?,n?,?b)=
        +Gtrace(?a,?b)*4*g(n1,n2)
        +rat(-2*ep,1)*Gtrace(?a,n1,n2,?b);
    id once ifmatch->retry Gtrace(?a,n?,n1?,n2?,n3?,n?,?b)=
        -2*Gtrace(?a,n3,n2,n1,?b)
        +rat(2*ep,1)*Gtrace(?a,n1,n2,n3,?b);
    id once ifmatch->retry Gtrace(?a,n?,n1?,n2?,n3?,n4?,n?,?b) =
        +2*Gtrace(?a,n4,n1,n2,n3,?b)
        +2*Gtrace(?a,n3,n2,n1,n4,?b)
        +rat(-2*ep,1)*Gtrace(?a,n1,n2,n3,n4,?b);
    id once ifmatch->retry Gtrace(?a,n?,n1?,...,n5?,n?,?b) =
        +2*Gtrace(?a,n5,n1,...,n4,?b)
        -2*Gtrace(?a,n4,n1,...,n3,n5,?b)
        -2*Gtrace(?a,n3,n2,n1,n4,n5,?b)
        +rat(2*ep,1)*Gtrace(?a,n1,...,n5,?b);
    id once ifmatch->retry Gtrace(?a,n?,n1?,...,n6?,n?,?b) =
        +2*Gtrace(?a,n6,n1,...,n5,?b)
        -2*Gtrace(?a,n5,n1,...,n4,n6,?b)
        +2*Gtrace(?a,n4,n1,...,n3,n5,n6,?b)
        +2*Gtrace(?a,n3,n2,n1,n4,n5,n6,?b)
        +rat(-2*ep,1)*Gtrace(?a,n1,...,n6,?b);
    id once ifmatch->retry Gtrace(?a,n?,n1?,...,n7?,n?,?b) =
        +2*Gtrace(?a,n7,n1,...,n6,?b)
        -2*Gtrace(?a,n6,n1,...,n5,n7,?b)
        +2*Gtrace(?a,n5,n1,...,n4,n6,n7,?b)
        -2*Gtrace(?a,n4,n1,...,n3,n5,n6,n7,?b)
        -2*Gtrace(?a,n3,n2,n1,n4,n5,n6,n7,?b)
        +rat(2*ep,1)*Gtrace(?a,n1,...,n7,?b);
* the generic fall-back case:
    id once ifmatch->retry Gtrace(?a,n?,?b,n1?,n?,?c) =
        +2*Gtrace(?a,n1,?b,?c) - Gtrace(?a,n,?b,n,n1,?c);
#endprocedure
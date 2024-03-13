Base.ceil(T::Type, x::Irrational) = Base.ceil(T, AbstractFloat(x))

function lower_cf(q::Number; maxitr = Inf)
  a=Int[]
  i=0
  while q != Inf && i < maxitr
    push!(a, q ÷ 1)
    q = 1/(q % 1)
    i=i+1
  end
  a
end

function upper_cf(q::Number; maxitr = Inf)
  a=Int[]
  i=0
  while q != 1 && i < maxitr
    push!(a, ceil(Int, q))
    q = 1/(1 - q%1)
    i=i+1
  end
  a
end

function farey!(a::Vector)
  for i in length(a):-1:2
    insert!(a,i,a[i] + a[i-1])
  end
  a
end

#=
function markov_extended_helper(p_1, q_1, p_2)
  p_3 = (3*p_1*p_2 + isqrt(9*p_1^2*p_2^2 - 4(p_1^2+p_2^2))) ÷ 2
  q_2 = (p_2 * q_1 -3p_3) ÷ p_1
  q_3 = (p_1*p_2*q_1*q_2 -q_2^2*p_1^2 -p_1*q_1-p_2*q_2) ÷ p_3

  q_2 = rem(q_2, p_2, RoundDown); q_2 = min(q_2,abs(-q_2))
  q_3 = rem(q_3, p_3, RoundDown); q_3 = min(q_3,abs(-q_3))

  q_2, p_3, q_3
end

function markov_extended(l,m,r)
  left = markov_extended_helper(l[1],l[2],m[1])
  ln = (l,[left[2],left[3]],m)

  right = markov_extended_helper(m[1],m[2],r[1])
  rn = (m,[right[2],right[3]],r)
  return (ln,rn)
end
=#


function neg_shear(v)
  -Rational[1-v[1]*v[2] v[1]^2;-v[2]^2 1+v[1]*v[2]]
end

function markov_extended(t)
  [[t[1], neg_shear(t[2])^(-1)*t[3], t[2]],
   [t[2], neg_shear(t[2])     *t[1], t[3]]]
end


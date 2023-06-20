function lower_cf(q::Rational)
  a=Int[]
  while q != 1//0
    push!(a, q ÷ 1)
    q = 1//(q % 1)
  end
  a
end

function upper_cf(q::Rational)
  a=Int[]
  while q != 1//1
    push!(a, ceil(Int, q))
    q = 1//(1-q % 1)
  end
  a
end

function fahry!(a::Vector)
  for i in length(a):-1:2
    insert!(a,i,a[i] .+ a[i-1])
  end
  a
end

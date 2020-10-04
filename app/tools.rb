module Clifford
  class Tools

    def self.random_number(bits)
      OpenSSL::BN::rand(bits).to_i
    end

    def self.random_prime(bits)
      # recall that OpenSLL library does not generate primes with less than
      # 16 bits
      if bits < 16
        Prime.first(100).select{|prime| prime.bit_length == bits}.sample
      else
        OpenSSL::BN::generate_prime(bits).to_i
      end
    end

    def self.next_prime(n)
      while true
        break if OpenSSL::BN.new(n).prime?
        n += 1
      end
       n
    end

    def self.generate_random_input(b)
      Array.new(8){ random_number(b) }
    end

    def self.generate_random_multivector_mod(b,q)
      data = generate_random_input(b)
      Clifford::Multivector3DMod.new data, q
    end

    def self.number_to_random_multivector_mod(n,b,q,g)
      e0,e1,e2,e3,e13,e23,e123 = Array.new(7){random_number(b)}
      e12 = -e0 - e1 + e2 + e3 - e13 + e23 + e123 + n
      data = [e0,e1,e2,e3,e12,e13,e23,e123]
      m_ = Clifford::Multivector3DMod.new data, q
      m = m_.scalar_mul(g)
      m
    end

    def self.exp_packing(n,b,q,g)
      n_ = random_number(g.bit_length - 2) * g + n
      input = (0..7).map{|i| (n_ / 2**(b*i)) % 2**b}
      m = Clifford::Multivector3DMod.new input, q
    end

    def self.exp_unpacking(m,b,q,g)
      n_ = 0; m.data.reverse.map{|d| n_ = n_ * 2**b + d}; n_
      n = n_ % g
    end

    def self.multivector_to_number(m,b,q,g)
      m_ = m.scalar_div(g)
      e0,e1,e2,e3,e12,e13,e23,e123 = m_.data
      n = e0 + e1 - e2 - e3 + e12 + e13 - e23 - e123
      n
    end

    def self.mod_inverse(num, mod)
      g, a, b = extended_gcd(num, mod)
      unless g == 1
        raise ZeroDivisionError.new("#{num} has no inverse modulo #{mod}")
      end
      a % mod
    end

    def self.extended_gcd(x, y)
      if x < 0
        g, a, b = extended_gcd(-x, y)
        return [g, -a, b]
      end
      if y < 0
        g, a, b = extended_gcd(x, -y)
        return [g, a, -b]
      end
      r0, r1 = x, y
      a0 = b1 = 1
      a1 = b0 = 0
      until r1.zero?
        q = r0 / r1
        r0, r1 = r1, r0 - q*r1
        a0, a1 = a1, a0 - q*a1
        b0, b1 = b1, b0 - q*b1
      end
      [r0, a0, b0]
    end

    def self.eea(p,m_10)
       a = [p, m_10]
       y = [0, 1]
       i = 1

       while a[i] > Integer.sqrt(p/2)
         q = a[i-1] / a[i]
         a << (a[i-1] - q * a[i])
         y << (y[i-1] + q * y[i])
         i += 1
       end

       Rational(a[i],y[i])
     end

     def self.modular_pow( base, power, mod )
      res = 1
      while power > 0
        res = (res * base) % mod if power & 1 == 1
        base = base ** 2 % mod
        power >>= 1
      end
      res
    end

  end
end

classdef tdist < pdfclass
  properties
    nu;
    mu;
    Sigma;

    MC = 1e7;
  end

  methods
    function this = tdist(mu, Sigma, nu)
      this.mu = mu;
      this.Sigma = Sigma;
      this.nu = nu;
    end

    function x = rand(this, n, m)
      if nargin < 3
        m = 1;
      end
      if nargin < 2
        n = 1;
      end
      y = sqrtm(this.Sigma) * randn(n, m);
      w = gamrnd(this.nu/2, 2/this.nu, n, m);
      x = bsxfun(@plus, this.mu, y./sqrt(w));
    end

    function Sigma = var(this)
      Sigma = this.nu/(this.nu-2)*this.Sigma;
    end

    function p = pdf(this, x)
      Delta = bsxfun(@minus, x, this.mu);
      Delta2 = sum(Delta.*(pinv(this.Sigma)*Delta), 1);
      p = exp(gammaln(this.nu/2+1)-gammaln(this.nu/2))/(this.nu*pi)^(size(this.mu, 1)/2)/...
          sqrt(det(this.Sigma))*(1+Delta2/this.nu).^-(this.nu/2+1);
    end

    function I = ia(this)
      x = rand(this, 1, this.MC);
      Delta = bsxfun(@minus, x, this.mu);
      Delta2 = sum(Delta.*(pinv(this.Sigma)*Delta), 1);

      I = (-(this.nu+2)/this.nu)^2 * mean(((pinv(this.Sigma)*x)./(1+Delta2/this.nu)).^2);
    end
  end
end

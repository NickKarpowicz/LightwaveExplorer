deviceFunction constexpr static inline uint32_t ufactorial(const uint32_t x){
    uint32_t f = 2u;
    if (x < 2) return 1u;
    for(uint32_t i = 3u; i<=x; i++){
        f *= i;
    }
    return f;
}

deviceFunction static inline deviceFP laguerre_prefactor(uint32_t p, uint32_t l){
    return deviceFPLib::sqrt(static_cast<deviceFP>((2 * ufactorial(p))/vPi<deviceFP>() * ufactorial(p + l)));
}

deviceFunction static inline deviceFP generalized_laguerre(const deviceFP x, const uint8_t alpha, const uint8_t n){
    if(n==0u) return 1.0f;
    deviceFP Lminus = static_cast<deviceFP>(1u + alpha) - x;
    if(n==1u) return Lminus;
    deviceFP Lminusminus = 1.0f;
    deviceFP Lk = {};
    const uint8_t alpha_plus_1 = alpha+1;

    for(uint8_t k = 1; k < n; k++){
        Lk = ((static_cast<deviceFP>(2u * k + alpha_plus_1) - x) * Lminus - static_cast<deviceFP>(k + alpha) * Lminusminus)/static_cast<deviceFP>(k + 1u);
        Lminusminus = Lminus;
        Lminus = Lk;
    }
    return Lk;
}

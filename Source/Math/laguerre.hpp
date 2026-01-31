deviceFunction static inline deviceFP laguerre_prefactor(const uint32_t p, const uint32_t l){
    deviceFP product = vPi<deviceFP>();
    for (unsigned int k = 1; k <= l; ++k) {
        product *= static_cast<deviceFP>(p + k);
    }
    return deviceFPLib::sqrt(2.0f/product);
}

deviceFunction static inline deviceFP generalized_laguerre(const deviceFP x, const uint32_t alpha, const uint32_t n){
    switch(n){
        case 0u: return 1.0f;
        case 1u: return static_cast<deviceFP>(1u + alpha) - x;
        case 2u: return 0.5f * (x * x - static_cast<deviceFP>(2u * (alpha + 2u)) * x + static_cast<deviceFP>((alpha + 1u)*(alpha + 2u)));
    }
    deviceFP Lk = static_cast<deviceFP>(1u + alpha) - x;
    deviceFP Lminus = 1.0f;
    const uint32_t alpha_plus_1 = alpha+1u;
    for(uint32_t k = 1; k < n; k++){
        deviceFP Lnext = ((static_cast<deviceFP>(2u * k + alpha_plus_1) - x) * Lk - static_cast<deviceFP>(k + alpha) * Lminus)/static_cast<deviceFP>(k + 1u);
        Lminus = Lk;
        Lk = Lnext;
    }
    return Lk;
}

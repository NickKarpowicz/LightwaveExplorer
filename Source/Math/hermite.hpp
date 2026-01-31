deviceFunction static inline constexpr deviceFP hermite(const deviceFP x, const uint32_t n){
    switch(n){
        case 0u: return 1.0f;
        case 1u: return 2.0f * x;
        case 2u: return 4.0f * x * x - 2.0f;
        case 3u: return 8.0f * x * x * x - 12.0f * x;
        default: {
            const deviceFP two_x = 2.0f*x;
            deviceFP Hminus = 1.0f;
            deviceFP H = two_x;
            for(uint32_t k = 1; k<n; k++){
                deviceFP Hnext = two_x * H - static_cast<deviceFP>(2 * k) * Hminus;
                Hminus = H;
                H = Hnext;
            }
            return H;
        }
    }
}

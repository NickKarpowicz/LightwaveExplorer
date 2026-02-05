deviceFunction static inline constexpr deviceFP hermite(const deviceFP x, const uint32_t n){
    switch(n){
        case 0u: return 1.0f;
        case 1u: return x;
        default: {
            deviceFP Hminus = 1.0f;
            deviceFP H = x;
            deviceFP normalization_factor = 1.0f;
            for(uint32_t k = 1; k<n; k++){
                deviceFP Hnext = x * H - static_cast<deviceFP>(k) * Hminus;
                //sqrt outside of loop would be faster but overflows too fast
                normalization_factor *= deviceFPLib::sqrt(static_cast<deviceFP>(k+1u));
                Hminus = H;
                H = Hnext;
            }
            return H/normalization_factor;
        }
    }
}

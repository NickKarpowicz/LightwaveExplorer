deviceFunction static inline constexpr deviceFP hermite(const deviceFP x, const unsigned int n){
    if(n==0) return 1.0f;
    const deviceFP two_x = 2.0f*x;
    if(n==1) return two_x;

    deviceFP Hminusminus = 1.0f;
    deviceFP Hminus = two_x;
    deviceFP H = {};

    for(unsigned int k = 1; k<n; k++){
        H = two_x * Hminus - static_cast<deviceFP>(2 * k) * Hminusminus;
        Hminusminus = Hminus;
        Hminus = H;
    }
    return H;
}

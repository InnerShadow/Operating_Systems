#include <iostream>
#include <math.h>

int nb_dig(const int n, const int d){
    int dig = 0;
    for(int k = 1; k <= n; ++k){
        int k_sq = std::pow(k, 2);

        while(k_sq){
            int current_dig = k_sq % 10;
            dig += (current_dig == d);
            k_sq /= 10;
        }
    }
    return dig;
}

int main(void){

    int n = 10;
    int d = 1;

    std::cout << (nb_dig(n, d)) << "\n";

    return 0;
}


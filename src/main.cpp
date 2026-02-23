#include "fp.h"
#include <iostream>

// demo
int main() {
    mpz_class p("7919");  // prime

    Fp a(42, p);
    Fp b(1337, p);

    std::cout << "p  = " << p << "\n";
    std::cout << "a  = " << a.val << "\n";
    std::cout << "b  = " << b.val << "\n";
    std::cout << "a+b = " << (a + b).val << "\n";
    std::cout << "a*b = " << (a * b).val << "\n";
    std::cout << "a/b = " << (a / b).val << "\n";
    std::cout << "a^(p-1) = " << pow(a, p - 1).val << "  (Fermat)\n";

    return 0;
}

#include <bits/stdc++.h>
using namespace std;

void solve(){
    double EPS = 1e-9;
    
    for (double x = -4; x < 4.1; x += 0.5) {
        double sum = 0;
        double term = 1;
        int n = 0;
        
        while (fabs(term) >= EPS) {
            sum += term;
            n += 1;
            term = term * (-x * x) / ((2.0 * n - 1.0) * (2.0 * n));
        }
        
        double exact = cos(x);
        double diff = fabs(sum - exact);

        cout << "x: " << x << " | sum: " << sum 
             << " | exact: " << exact 
             << " | diff: " << diff << '\n';
    }
}

int main(){
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    cout.tie(0);
    solve();
    return 0;
}

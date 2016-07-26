#include <iostream>

using namespace std;

int main() {
    int volumeThickness;
    cin >> volumeThickness;
    int coverThickness;
    cin >> coverThickness;
    int start;
    cin >> start;
    int end;
    cin >> end;
    if (end > start) {
        cout << 2*coverThickness + (end - start - 1)*(volumeThickness + 2*coverThickness) << endl;
    } else if (start == end) {
        cout << volumeThickness << endl;
    } else {
        cout << volumeThickness + (start - end)*(volumeThickness + 2*coverThickness) << endl;
    }
    return 0;
}

#include <nsga2/global.h>
#include <vector>
#include <iostream>

using namespace std;


int main(int argc, char *argv[]) {
    std::vector<int> nbits(3,5);
    nsga2::individual ind(3,3,2,nbits,2);

    cout << "Individual: " << ind << endl;
    
    return 0;
}

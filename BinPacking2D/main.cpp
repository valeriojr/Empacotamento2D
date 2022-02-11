#include <iostream>

#include "BRKGA.h"
#include "MTRand.h"
#include "3rdParty/RectangleBinPack/SkylineBinPack.h"
#include "3rdParty/RectangleBinPack/ShelfBinPack.h"

#include <RectangleBinPack/Rect.h>
#include <RectangleBinPack/GuillotineBinPack.h>
#include <RectangleBinPack/MaxRectsBinPack.h>


template<class T, typename H>
class BP2Decoder {
private:
    const double width;
    std::vector<rbp::Rect>& boxes;
    T* packer;
    H heuristic;
public:
    BP2Decoder(const double width, std::vector<rbp::Rect>& boxes, const H& heuristic) : boxes(boxes), width(width),
                                                                                        heuristic(heuristic) {
        packer = new T();
    }

    ~BP2Decoder() {
        delete packer;
    }

    double decode(const std::vector<double>& chromosome) const {
        packer->Init(width, INT_MAX);

        std::vector<std::pair<double, int>> ordering(chromosome.size());
        for(int i = 0;i < chromosome.size();i++)
            ordering[i].second = i;

        std::sort(ordering.begin(), ordering.end());

        int height = 0;
        for (auto & box : ordering) {
            rbp::Rect placed = packer->Insert(boxes[box.second].width, boxes[box.second].height, heuristic);
            height = std::max(height, placed.y + placed.height);
        }

        return -height;
    }
};

typedef BP2Decoder<rbp::MaxRectsBinPack, rbp::MaxRectsBinPack::FreeRectChoiceHeuristic> MaxRectsDecoder;
typedef BP2Decoder<rbp::GuillotineBinPack, rbp::GuillotineBinPack::FreeRectChoiceHeuristic> GuillotineDecoder;
typedef BP2Decoder<rbp::ShelfBinPack, rbp::ShelfBinPack::ShelfChoiceHeuristic> ShelfDecoder;
typedef BP2Decoder<rbp::SkylineBinPack, rbp::SkylineBinPack::LevelChoiceHeuristic> SkyLineDecoder;


int main() {
    unsigned n = 10; // chromosome size (number of boxes)
    unsigned p = 1000; // population size
    const double pe = 0.1;
    const double pm = 0.1;
    const double rhoe = 0.7;
    const unsigned K = 1;
    double width;

    std::cout << "Enter the strip's width" << std::endl;
    std::cin >> width;
    std::cout << "Enter the number of boxes" << std::endl;
    std::cin >> n;
    std::vector<rbp::Rect> boxes(n);
    for (int i = 0;i < n;i++) {
        std::cout << "Enter box #" << i << " dimensions" << std::endl;
        std::cin >> boxes[i].width >> boxes[i].height;
    }

    MaxRectsDecoder decoder(width, boxes, rbp::MaxRectsBinPack::RectBottomLeftRule);
    MTRand rng;

    BRKGA<MaxRectsDecoder, MTRand> brkga(n, p, pe, pm, rhoe, decoder, rng);

    unsigned generation = 0;        // current generation
    const unsigned X_INTERVAL = 100;    // exchange best individuals at every 100 generations
    const unsigned X_NUMBER = 2;    // exchange top 2 best
    const unsigned MAX_GENS = 1000;    // run for 1000 gens
    std::cout << "Running for " << MAX_GENS << " generations..." << std::endl;
    do {
        // Salvar solução inicial
        brkga.evolve();    // evolve the population for one generation

        if ((++generation) % X_INTERVAL == 0) {
            brkga.exchangeElite(X_NUMBER);    // exchange top individuals
        }
    } while (generation < MAX_GENS);

    // print the fitness of the top 10 individuals of each population:
    std::cout << "Fitness of the top 10 individuals of each population:" << std::endl;
    const unsigned bound = std::min(p, unsigned(10));    // makes sure we have 10 individuals
    for (unsigned i = 0;i < K;++i) {
        std::cout << "Population #" << i << ":" << std::endl;
        for (unsigned j = 0;j < bound;++j) {
            std::cout << "\t" << j << ") "
                      << brkga.getPopulation(i).getFitness(j) << std::endl;
        }
    }

    std::cout << "Best solution found has objective value = "
              << brkga.getBestFitness() << std::endl;

    return 0;
}

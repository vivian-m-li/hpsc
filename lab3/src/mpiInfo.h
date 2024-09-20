
//  ==
//  ||
//  ||   C L A S S:    m p i I n f o
//  ||
//  ==

class mpiInfo {
 public:
  int myPE;
  int numPE;
  int nRealx, nRealy;
  int nx, ny;
  int nPEx, nPEy;
  int iPE, jPE;
  int iMin, iMax, jMin, jMax;  // The global i-j numbers on this processor
  int nei_n, nei_s, nei_e, nei_w;
  int nei_ne, nei_nw, nei_se, nei_sw;

  int countx, county;

  double *phiL, *phiR;
  double *phiT, *phiB;

  double *phiSend_n, *phiSend_s;
  double *phiSend_e, *phiSend_w;
  double *phiRecv_n, *phiRecv_s, *phiRecv_e, *phiRecv_w;

  MPI_Status status;
  int err;
  int tag;
  MPI_Request request;

  //  -
  //  |
  //  |   GridDecomposition: Set up PE numbering system in figure below and
  //  |                      establish communication arrays.
  //  |
  //  |                      nPEx -- number of PEs in the x-direction
  //  |                      nPEy -- number of PEs in the y-direction
  //  |                      numPE = total number of PEs
  //  |
  //  |                       +-------+-------+ . . .   +-------+
  //  |                       |       |       |         |       |
  //  |                       |       |       |         | numPE |
  //  |                       |       |       |         |       |
  //  |                       +-------+-------+ . . .   +-------+
  //  |                       .       .       .         .       .
  //  |                       .       .       .         .       .
  //  |                       .       .       .         .       .
  //  |                       +-------+-------+ . . .   +-------+
  //  |                       |       |       |         |       |
  //  |                       | nPEx  | nPEx+1|         |       |
  //  |                       |       |       |         |       |
  //  |                       +-------+-------+ . . .   +-------+
  //  |                       |       |       |         |       |
  //  |                       |   0   |   1   |         | nPEx-1|
  //  |                       |       |       |         |       |
  //  |                       +-------+-------+ . . .   +-------+
  //  |
  //  |
  //  -

  void GridDecomposition(int _nPEx, int _nPEy, int nCellx, int nCelly) {
    nRealx = nCellx + 1;
    nRealy = nCelly + 1;

    // Store and check incoming processor counts

    nPEx = _nPEx;
    nPEy = _nPEy;

    if (nPEx * nPEy != numPE) {
      if (myPE == 0)
        cout << "Fatal Error:  Number of PEs in x-y directions do not add up "
                "to numPE"
             << endl;
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Finalize();
      exit(0);
    }

    // Get the i-j location of this processor, given its number.  See figure
    // above:

    jPE = int(myPE / nPEx);
    iPE = myPE - jPE * nPEx;

    // Set neighbor values

    nei_n = nei_s = nei_e = nei_w = -1;

    if (iPE > 0) {
      nei_w = myPE - 1;
    }
    if (jPE > 0) {
      nei_s = myPE - nPEx;
    }
    if (iPE < nPEx - 1) {
      nei_e = myPE + 1;
    }
    if (jPE < nPEy - 1) {
      nei_n = myPE + nPEx;
    }

    nei_nw = nei_sw = nei_ne = nei_se = -1;

    if (iPE > 0 && jPE > 0) nei_sw = /* TO-DO in Lab */;
    if (iPE < nPEx - 1 && jPE > 0) nei_se = /* TO-DO in Lab */;
    if (iPE > 0 && jPE < nPEy - 1) nei_nw = /* TO-DO in Lab */;
    if (iPE < nPEx - 1 && jPE < nPEy - 1) nei_ne = /* TO-DO in Lab */;

    // Acquire memory for the communication between adjacent processors:
    countx = nRealx + 2;
    county = nRealy + 2;

    phiL = new double[county];
    phiR = new double[county];
    phiT = new double[countx];
    phiB = new double[countx];

    phiSend_n = new double[countx];
    phiSend_s = new double[countx];
    phiSend_e = new double[county];
    phiSend_w = new double[county];
    phiRecv_n = new double[countx];
    phiRecv_s = new double[countx];
    phiRecv_e = new double[county];
    phiRecv_w = new double[county];

    tag = 0;
  }

  //  ==
  //  ||
  //  ||  ParticlesExchange
  //  ||
  //  ||  Exchange particles between processors
  //  ||
  //  ||
  //  ==

  void ParticleExchange(VI &ptcl_send_list, VI &ptcl_send_PE, particles &PTCL) {
    MPI_Status status;
    MPI_Request request;

    // (1) Get the max number particles to be sent by any particular processor,
    // and make sure all processors  know that number.

    int numToSend = ptcl_send_list.size();
    int maxToSend;

    MPI_Iallreduce(/* TO-DO in Lab */, /* TO-DO in Lab */, 1, MPI_INT,
                   /* TO-DO in Lab */, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);

    // (2) Allocate contributions to the upcoming Gather operation.  Here, "C"
    // for "Contribution" to be Gathered

    int *Cptcl_PE;
    Cptcl_PE = new int[/* TO-DO in Lab */];  // Particles' destination PEs
    double *Cptcl_x;
    Cptcl_x = new double[/* TO-DO in Lab */];
    double *Cptcl_y;
    Cptcl_y = new double[/* TO-DO in Lab */];
    double *Cptcl_vx;
    Cptcl_vx = new double[/* TO-DO in Lab */];
    double *Cptcl_vy;
    Cptcl_vy = new double[/* TO-DO in Lab */];

    // (3) Populate contributions on all processors for the upcoming Gather
    // operation

    for (int i = 0; i < maxToSend; ++i) {
      Cptcl_PE[i] = -1;
      Cptcl_x[i] = 0.;
      Cptcl_y[i] = 0.;
      Cptcl_vx[i] = 0.;
      Cptcl_vy[i] = 0.;
    }

    // (4) Populate with all the particles on this PE.  Note that some/most
    // processors will have left-over space in the C* arrays.

    for (int i = 0; i < ptcl_send_list.size(); ++i) {
      int id = ptcl_send_list[/* TO-DO in Lab */];
      Cptcl_PE[i] = ptcl_send_PE[/* TO-DO in Lab */];
      Cptcl_x[i] = PTCL.x[/* TO-DO in Lab */];
      Cptcl_y[i] = PTCL.y[/* TO-DO in Lab */];
      Cptcl_vx[i] = PTCL.vx[/* TO-DO in Lab */];
      Cptcl_vy[i] = PTCL.vy[/* TO-DO in Lab */];
    }

    // (5) Allocate and initialize the arrays for upcoming Gather operation to
    // PE0.  The sizeOfGather takes
    //     into account the number of processors, like this figure:
    //
    //     |<----------------------------- sizeOfGather
    //     ------------------------------>| | | | |
    //     |<- maxToSend    ->|<- maxToSend    ->|<- maxToSend    ->|<-
    //     maxToSend    ->|
    //     +------------------+------------------+------------------+------------------+
    //             PE0               PE1                PE2               PE3

    int sizeOfGather = /* TO-DO in Lab */;

    int *Gptcl_PE;
    Gptcl_PE = new int[/* TO-DO in Lab */];
    double *Gptcl_x;
    Gptcl_x = new double[/* TO-DO in Lab */];
    double *Gptcl_y;
    Gptcl_y = new double[/* TO-DO in Lab */];
    double *Gptcl_vx;
    Gptcl_vx = new double[/* TO-DO in Lab */];
    double *Gptcl_vy;
    Gptcl_vy = new double[/* TO-DO in Lab */];

    for (int i = 0; i < sizeOfGather; ++i) {
      Gptcl_PE[i] = -1;
      Gptcl_x[i] = 0.;
      Gptcl_y[i] = 0.;
      Gptcl_vx[i] = 0.;
      Gptcl_vy[i] = 0.;
    }

    // (6)  Gather "Contributions" ("C" arrays) from all PEs onto all PEs into
    // these bigger arrays so all PE will know what particles
    //      need to go where.

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Iallgather(/*  TO-DO in Lab */, &request);
    MPI_Wait(&request, &status);
    MPI_Iallgather(/*  TO-DO in Lab */, &request);
    MPI_Wait(&request, &status);
    MPI_Iallgather(/*  TO-DO in Lab */, &request);
    MPI_Wait(&request, &status);
    MPI_Iallgather(/*  TO-DO in Lab */, &request);
    MPI_Wait(&request, &status);
    MPI_Iallgather(/*  TO-DO in Lab */, &request);
    MPI_Wait(&request, &status);

    MPI_Barrier(MPI_COMM_WORLD);

    // (7) Put in vector form so they can be added to PTCL.  These arrays are
    // 1-based.

    int Np = 0;
    for (int i = 0; i < /*  TO-DO in Lab */; ++i)
      if (Gptcl_PE[i] == myPE) ++Np;

    VD std_add_x;
    std_add_x.resize(Np + 1);
    VD std_add_y;
    std_add_y.resize(Np + 1);
    VD std_add_vx;
    std_add_vx.resize(Np + 1);
    VD std_add_vy;
    std_add_vy.resize(Np + 1);

    int count = 1;
    for (int i = 0; i < sizeOfGather; ++i)
      if (Gptcl_PE[i] == /* TO-DO in Lab */) {
        std_add_x[/* TO-DO in Lab */] = Gptcl_x[i];
        std_add_y[/* TO-DO in Lab */] = Gptcl_y[i];
        std_add_vx[/* TO-DO in Lab */] = Gptcl_vx[i];
        std_add_vy[/* TO-DO in Lab */] = Gptcl_vy[i];
        /* TO-DO in Lab */;
      }

    PTCL.add(/* TO-DO in Lab */);

    // (8) Free up memory

    if (maxToSend > 0) {
      delete[] Cptcl_PE;
      delete[] Cptcl_x;
      delete[] Cptcl_y;
      delete[] Cptcl_vx;
      delete[] Cptcl_vy;
    }
    if (sizeOfGather > 0) {
      delete[] Gptcl_PE;
      delete[] Gptcl_x;
      delete[] Gptcl_y;
      delete[] Gptcl_vx;
      delete[] Gptcl_vy;
    }
  }

  int pid(int i, int j) { return (i + 1) + (j) * (nRealx + 2); }
};

/*
 * cache_controller.c
 *
 *  Created on: Nov 29, 2022
 *      Author: madhu,nithin
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>

#define complex _Complex
#define i _Imaginary_I
#define CMPLX(x, y) ((double complex)((double)(x) + _Imaginary_I * (double)(y)))

#define MAX_LINES 65536
#define MAX_N_WAY 16

#define CACHE_NOT_HIT -1
#define ALL_VALID     -1

//write strategy
#define WB  0
#define WTA  1
#define WTNA  2

uint32_t old_line;
uint32_t g_line;
uint32_t g_tag;
uint32_t cache_size = (256*1024);
uint32_t addrBits = 32;
uint32_t cacheLine, l, b, t, block_size, total_lines, N_WAY;
uint8_t dataBits = 32;
uint8_t dataBytes = 4;

bool cacheValid[MAX_LINES][MAX_N_WAY];
bool cacheDirty[MAX_LINES][MAX_N_WAY];
uint8_t cacheLRU[MAX_LINES][MAX_N_WAY];
uint32_t cacheTag[MAX_LINES][MAX_N_WAY];


//counters for read
uint32_t rmc, rbc, rbhc, rbmc, rbrdc, rbrc;

//counters for write
uint32_t wmc, wbc, wtc, wbhc, wbmc, wbrdc, wbrc;

uint32_t flush_count;

void resetCache()
{
    rmc = 0;
    rbc = 0;
    rbhc = 0;
    rbmc = 0;
    rbrdc = 0;
    rbrc = 0;

    wmc = 0;
    wbc = 0;
    wtc = 0;
    wbhc = 0;
    wbmc = 0;
    wbrdc = 0;
    wbrc = 0;

    flush_count = 0;
    for(int i = 0; i < MAX_LINES; i++)
    {
        for(int j=0; j < MAX_N_WAY; j++)
        {
            cacheDirty[i][j] = false;
            cacheValid[i][j] = false;
            cacheTag[i][j] = 0;
            cacheLRU[i][j] = j;
        }
    }
}

void flush_cache()
{
    for(int i = 0; i < MAX_LINES; i++)
    {
        for(int j = 0; j < MAX_N_WAY; j++)
        {
            if((cacheValid[i][j] == true) && (cacheDirty[i][j] == true))
            {
                flush_count++;
                cacheDirty[i][j] = false;
            }
        }
    }

}

int ifCacheHit(uint32_t line, uint32_t tag)
{
    int i=0;
    while(i<N_WAY)
    {
        if(cacheTag[line][i] == tag)
        {
            if(cacheValid[line][i] == true)
                return i;
        }
        i++;
    }
    return CACHE_NOT_HIT;
}

int allBlocksValid(uint32_t line)
{
    int i = 0;
    while(i<N_WAY)
    {
        if(cacheValid[line][i] == 0)
            return i;
        i++;
    }
    return ALL_VALID;
}

uint32_t fnLine(uint8_t *address)
{
    return (((uint32_t)address<<t)>>(t+b));
}

uint32_t fnTag(uint8_t *address)
{
    return ((uint32_t)address>>(l+b));
}

void update_LRU(uint32_t line, uint8_t hit_LRU)
{
    int i = 0;
    for(i = 0; i<N_WAY; i++)
    {
        if(cacheLRU[line][i] < hit_LRU)
            cacheLRU[line][i]++;
    }
    cacheLRU[line][hit_LRU] = 0;
}

int oldestLRU(uint32_t line)
{
    int i;
    for(i = 0; i<N_WAY; i++)
    {
        if(cacheLRU[line][i] == (N_WAY-1))
        {
            return i;
        }
    }
}

void read_Block(uint32_t line, uint32_t tag)
{
    rbc++;
    int hit_LRU, rd_index = 0;
    rd_index = ifCacheHit(line,tag);
    if(rd_index != CACHE_NOT_HIT)
    {
        rbhc++;
    }
    else
    {
        rbmc++;
        rd_index = allBlocksValid(line);
        if(rd_index == ALL_VALID)
        {
            rd_index = oldestLRU(line);
            if(cacheDirty[line][rd_index] == true)
                {
                  rbrdc++;
                  cacheDirty[line][rd_index] = false;
                }
            rbrc++;
        }
      }
    if(rd_index != ALL_VALID)
    {
    			cacheValid[line][rd_index] = true;
    	        cacheDirty[line][rd_index] = false;
    	        cacheTag[line][rd_index] = tag;
    	        hit_LRU = cacheLRU[line][rd_index];
    	        update_LRU(line, hit_LRU);
    }
}


void read_Memory(uint8_t *address, uint8_t byteCount)
{
    int i = 0;
    old_line = -9999;
    int line=0, tag=0;
    rmc++;
    for(i = 0; i<byteCount; i++)
    {
        line = fnLine(address);
        tag = fnTag(address);

        if(line != old_line)
        {
            read_Block(line, tag);
            old_line = line;
        }
    address++;
    }
}

void write_Block(uint32_t line, uint32_t tag, uint8_t WS)
{
    int wr_index = 0, hit_LRU;
    wbc++;
    wr_index = ifCacheHit(line,tag);
    if(wr_index != CACHE_NOT_HIT)
    {
        wbhc++;
    }
    else
    {
        wbmc++;
        wr_index = allBlocksValid(line);
        if((wr_index == ALL_VALID ) && !(WS == WTNA))
        {
            wr_index = oldestLRU(line);
            if(cacheDirty[line][wr_index] == true)
                {
                  wbrdc++;
                  cacheDirty[line][wr_index]=false;
                }
            wbrc++;
        }
    }
    if(wr_index!=ALL_VALID)
    {
        hit_LRU = cacheLRU[line][wr_index];
        update_LRU(line, hit_LRU);
        if(WS == WB)
        {
            cacheDirty[line][wr_index] = true;
        }
        cacheValid[line][wr_index] = true;
        cacheTag[line][wr_index] = tag;
    }
}


void write_Memory(uint8_t* address, uint8_t byteCount, uint8_t WS)
{
    wmc++;

    int i = 0, line=0, tag=0;
    old_line = -9999;
    for(i = 0; i < byteCount; i++)
    {
         line = fnLine(address);
         tag = fnTag(address);

        if(line != old_line)
        {
            if(( WS == WTA) || (WS == WTNA))
                wtc++;
            write_Block(line, tag, WS);
            old_line = line;
        }
        address++;
    }
}

void Radix2FFT(complex data[], int nPoints, int nBits, uint8_t WS)
{
  // cooley-tukey radix-2, twiddle factor
  // adapted from Fortran code from Burrus, 1983
  #pragma warning (disable: 4270)
  int i, j, k, l;
  int nPoints1, nPoints2;
  complex cTemp, cTemp2;
  double dTheta, dDeltaCos,dDeltaSin, dCos, dSin;

  nPoints2 = nPoints;
  read_Memory(&nPoints, sizeof(nPoints));
  write_Memory(&nPoints2, sizeof(nPoints2), WS);

  write_Memory(&k, sizeof(k),WS);
  read_Memory(&nBits, sizeof(nBits));

  for (k = 1; k <= nBits; k++)
  {
    read_Memory(&k, sizeof(k));
    read_Memory(&nBits, sizeof(nBits));

    nPoints1 = nPoints2;
    read_Memory(&nPoints2, sizeof(nPoints2));
    write_Memory(&nPoints1, sizeof(nPoints1), WS);

    nPoints2 /= 2;
    read_Memory(&nPoints2, sizeof(nPoints2));
    write_Memory(&nPoints2, sizeof(nPoints2), WS);


    // Compute differential angles
    double dTheta = 2 * 3.14159257 / nPoints1;
    read_Memory(&nPoints1, sizeof(nPoints1));
    write_Memory(&dTheta, sizeof(dTheta), WS);

    double dDeltaCos = cos(dTheta);
    read_Memory(&dTheta, sizeof(dTheta));
    write_Memory(&dDeltaCos, sizeof(dDeltaCos), WS);



    double dDeltaSin = sin(dTheta);
    read_Memory(&dTheta, sizeof(dTheta));
    write_Memory(&dDeltaSin, sizeof(dDeltaSin), WS);


    // Initialize angles
     double dCos = 1;
    write_Memory(&dCos, sizeof(dCos), WS);

    double dSin = 0;
    write_Memory(&dSin, sizeof(dSin), WS);

    // Perform in-place FFT
    write_Memory(&j, sizeof(j), WS);
    read_Memory(&nPoints2, sizeof(nPoints2));
    for (j = 0; j < nPoints2; j++)
    {
      read_Memory(&j, sizeof(j));
      read_Memory(&nPoints2, sizeof(nPoints2));

      i = j;
      read_Memory(&j, sizeof(j));
      write_Memory(&i, sizeof(i), WS);

      read_Memory(&i, sizeof(i));
      read_Memory(&nPoints, sizeof(nPoints));

      while (i < nPoints)
      {

        read_Memory(&i, sizeof(i));
        read_Memory(&nPoints, sizeof(nPoints));

        read_Memory(&i, sizeof(i));
        read_Memory(&nPoints2, sizeof(nPoints2));

        l = i + nPoints2;
        write_Memory(&l, sizeof(l), WS);

        cTemp = data[i] - data[l];
        read_Memory(&i, sizeof(i));
        read_Memory(&l, sizeof(l));
        read_Memory(&data[i], sizeof(data[i]));
        read_Memory(&data[l], sizeof(data[l]));
        write_Memory(&cTemp, sizeof(cTemp), WS);

        cTemp2 = data[i] + data[l];
        read_Memory(&i, sizeof(i));
        read_Memory(&l, sizeof(l));
        read_Memory(&data[i], sizeof(data[i]));
        read_Memory(&data[l], sizeof(data[l]));
        write_Memory(&cTemp2, sizeof(cTemp2), WS);

        data[i] = cTemp2;
        read_Memory(&cTemp2, sizeof(cTemp2));
        read_Memory(&i, sizeof(i));
        write_Memory(&data[i], sizeof(data[i]), WS);

        cTemp2 = CMPLX(dCos * creal(cTemp) + dSin * cimag(cTemp),dCos * cimag(cTemp) - dSin * creal(cTemp));
        read_Memory(&cTemp, sizeof(creal(cTemp)));
        read_Memory(&dCos, sizeof(dCos));
        read_Memory(&cTemp, sizeof(cimag(cTemp)));
        read_Memory(&dSin, sizeof(dSin));
        read_Memory(&cTemp, sizeof(creal(cTemp)));
        read_Memory(&dCos, sizeof(dCos));
        read_Memory(&cTemp, sizeof(cimag(cTemp)));
        read_Memory(&dSin, sizeof(dSin));
        write_Memory(&cTemp2, sizeof(cTemp2), WS);

        data[l] = cTemp2;
        read_Memory(&cTemp2, sizeof(cTemp2));
        read_Memory(&l, sizeof(l));
        write_Memory(&data[l], sizeof(data[l]),WS);

        i += nPoints1;
        read_Memory(&nPoints1, sizeof(nPoints1));
        read_Memory(&i, sizeof(i));
        write_Memory(&i, sizeof(i), WS);
        read_Memory(&i, sizeof(i));
        read_Memory(&nPoints, sizeof(nPoints));
      }

      double dTemp = dCos;
      read_Memory(&dCos, sizeof(dCos));
      write_Memory(&dTemp, sizeof(dTemp), WS);



      dCos = dCos * dDeltaCos - dSin * dDeltaSin;
      read_Memory(&dCos, sizeof(dCos));
      read_Memory(&dDeltaCos, sizeof(dDeltaCos));
      read_Memory(&dSin, sizeof(dSin));
      read_Memory(&dDeltaSin, sizeof(dDeltaSin));
      write_Memory(&dCos, sizeof(dCos), WS);

      dSin = dTemp * dDeltaSin + dSin * dDeltaCos;
      read_Memory(&dTemp, sizeof(dTemp));
      read_Memory(&dDeltaSin, sizeof(dDeltaSin));
      read_Memory(&dSin, sizeof(dSin));
      read_Memory(&dDeltaCos, sizeof(dDeltaCos));
      write_Memory(&dSin, sizeof(dSin), WS);

      write_Memory(&j, sizeof(j),WS);
    }

    write_Memory(&k, sizeof(k), WS);
  }

  //Convert Bit Reverse Order to Normal Ordering
    j = 0;
    write_Memory(&j, sizeof(j), WS);
    nPoints1 = nPoints - 1;
    read_Memory(&nPoints, sizeof(nPoints));
    write_Memory(&nPoints1, sizeof(nPoints1), WS);

    write_Memory(&i, sizeof(i), WS);
    read_Memory(&nPoints1, sizeof(nPoints1));
  for(i = 0; i < nPoints1; i++)
  {
       read_Memory(&i, sizeof(i));
       read_Memory(&nPoints1, sizeof(nPoints1));

       read_Memory(&i, sizeof(i));
       read_Memory(&j, sizeof(j));

       if (i < j) {
        cTemp = data[j];
        read_Memory(&j, sizeof(j));
        read_Memory(&data[j], sizeof(data[j]));
        write_Memory(&cTemp, sizeof(cTemp),WS);

        cTemp2 = data[i];
        read_Memory(&i, sizeof(i));
        read_Memory(&data[i], sizeof(data[i]));
        write_Memory(&cTemp2, sizeof(cTemp2), WS);

        data[i] = cTemp;
        read_Memory(&i, sizeof(i));
        read_Memory(&cTemp, sizeof(cTemp));
        write_Memory(&data[i], sizeof(data[i]), WS);

        data[j] = cTemp2;
        read_Memory(&j, sizeof(j));
        read_Memory(&cTemp2, sizeof(cTemp2));
        write_Memory(&data[j], sizeof(data[j]), WS);
        }
    k = nPoints / 2;
    read_Memory(&nPoints, sizeof(nPoints));
    write_Memory(&k, sizeof(k), WS);


    read_Memory(&j, sizeof(j));
    read_Memory(&k, sizeof(k));
    while (k <= j)
    {

        read_Memory(&k, sizeof(k));
        read_Memory(&j, sizeof(j));


    j -= k;
    read_Memory(&k, sizeof(k));
    read_Memory(&j, sizeof(j));
    write_Memory(&j, sizeof(j),WS);


    k /= 2;
    read_Memory(&k, sizeof(k));
    write_Memory(&k, sizeof(k), WS);
    read_Memory(&j, sizeof(j));
    }

    j += k;
    read_Memory(&k, sizeof(k));
    read_Memory(&j, sizeof(j));
    write_Memory(&j, sizeof(j), WS);

    write_Memory(&i, sizeof(i), WS);
    }
  #pragma warning(default: 4270)
}



int main(int argc, char* argv[])
{
  int i,j;
  double complex data[32768];

  FILE *fp;

  printf("FFT Test App\r\n");

  int points = 32768;
  // time domain: zero-offset real cosine function
  // freq domain: delta fn with zero phase
  #define cycles 1 // max cycles is points/2 for Nyquist
  for (i = 0; i < points; i++)
  {
      double temp = (cos(2.0*3.1416*(float)i/(float)(points)*cycles));
      data[i] = CMPLX(temp, 0.0);
  }

  int bits = ceil(log(points)/log(2));
  int BL, N, WS;
    printf("FFT Test App\r\n");
  fp = fopen("/Users/madhu/eclipse-workspace/6313_project/t8.csv", "w");

  fprintf(fp, "BL,N,WS,rmc, rbc, rbhc, rbmc, rbrdc, rbrc, wmc, wbc, wtc, wbhc, wbmc, wbrdc, wbrc,flushCount\n");

  for(N=0;N<=16;N++)
  {
    if(N==1 || N==2 || N==4 || N==8 || N==16)
    {
        for(BL = 1;BL<=8; BL++)
        {
        	if(BL==1 || BL==2 || BL==4 || BL==8)
        	{
            for(WS=0; WS<=2; WS++)
            {
            	N_WAY = N;
            	block_size = BL*dataBytes;
                b = log2(block_size);
                total_lines = (cache_size/(block_size*N));
                l = log2(total_lines);
            	t = dataBits-b-l;
                resetCache();

                Radix2FFT(data, points, bits, WS);
                flush_cache();
                fprintf(fp, "%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n",BL,N,WS,rmc, rbc, rbhc, rbmc, rbrdc, rbrc, wmc, wbc, wtc, wbhc, wbmc, wbrdc,wbrc,flush_count);
            }
            }
        }
    }
  }


  return 0;
}



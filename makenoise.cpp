//******************************************************************************
//    Copyright (c) 2010, Christopher James Huff
//    All rights reserved.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//******************************************************************************

// g++ -std=c++1y -O3 -Ilibs/include makenoise.cpp libs/lib/libFLAC++.la libs/lib/libFLAC.la -o makenoise && ./makenoise && open flacout.flac
// clang++ -std=c++1y -stdlib=libc++ -O3 makenoise.cpp -lFLAC -lFLAC++ -o makenoise && ./makenoise

#include <iostream>
#include <sys/time.h>

#include "soundbuffer.h"

using namespace std;

inline double GetRealSeconds() {
    timeval newTime;
    gettimeofday(&newTime, NULL);
    return newTime.tv_sec + newTime.tv_usec/1e6;
}

#if(0)
// Waveforms normally have a range of -1..1. Amplitude is determined by modulating
// with a decay function, use NoDecay() for constant amplitude.
struct SineWave {
    double phase, freq;// phase normalized from 0-1, frequency given in Hz and stored in radians/s
    SineWave(double f, double p): phase(p*pi2), freq(f*pi2) {}
    float operator()(float t) const {return sin(t*freq + phase);}
};

template<typename Fn>
struct NoDecay {
    Fn fn;
    float amplitude;
    NoDecay(const Fn & f, float a = 1.0): fn(f), amplitude(a) {}
    float operator()(float t) const {return fn(t)*amplitude;}
};
struct BoxcarDecay {
    float amplitude, duration;
    BoxcarDecay(float a, float dur): amplitude(a), duration(dur) {}
    float operator()(float t) const {return (t < duration)? amplitude : 0.0;}
};

template<typename Fn>
struct ExpDecay {
    Fn fn;
    float amplitude;
    float decay;
    ExpDecay(const Fn & f, float a, float d): fn(f), amplitude(a), decay(d) {}
    float operator()(float t) const {return fn(t)*amplitude*exp(-t*decay);}
};


// Sum sound function to sound buffer
template<typename Fn>
inline void AddSound(sndbfr_48k_t & sb, const Fn & fn) {
    for(size_t j = 0, nsamps = sb.size(); j < nsamps; ++j) {
        sb.data_l[j] += fn(sb.t(j));
        sb.data_r[j] = sb.data_l[j];
    }
}
template<typename Fn>
inline void AddSound(sndbfr_48k_t & sb, float startT, float endT, const Fn & fn) {
    size_t end = min((int)sb.size(), sb.index(endT));
    for(size_t j = sb.index(startT); j < end; ++j) {
        sb.data_l[j] += fn((float)j/FS);
        sb.data_r[j] = sb.data_l[j];
    }
}

void AddSine(sndbfr_48k_t & sb, float freq, float phase, float amp) {
    AddSound(sb, NoDecay<SineWave>(SineWave(freq, phase), amp));
}


template<typename F>
void Generate(sndbfr_48k_t & sb, int channel, float startT, float endT)
{
    size_t nsamps = sb.size();
    for(int j = sb.index(startT); j < nsamps; ++j) {
        sb[channel][j] += a*F(sb.t(j));
    }
}
#endif

// Add sine with exponential decay. Rough approximation of damped harmonic oscillator.
// Note that it won't be continuous if phase is not zero.
/*void AddDecaySine(sndbfr_48k_t & sb, sec_t startT, sec_t freq, sec_t phase, samp_t amp, sec_t decay)
{
    auto nsamps = sb.size();
    for(int j = sb.index(startT); j < nsamps; ++j) {
        sec_t t = sb.t(j) - startT;
        samp_t a = amp*exp(-t*decay);
        if(a < 1.0e-6f) break;
        sb[j] += a*sin(t*pi2*freq + phase);
    }
}

// Does an in-place low-pass, averaging n=taps samples into each output sample
// Leading samples are effectively averaged with zero samples.
void LowPass(sndbfr_48k_t & sb, int taps)
{
    auto nsamps = sb.size();
    for(int j = nsamps - 1; j > 0; --j) {
        for(int k = j - 1; k > 0 && k > (j - taps); --k)
            sb[j] += sb[k];
        sb[j] /= taps;
    }
}
void HighPass(sndbfr_48k_t & sb, int taps)
{
    sndbfr_48k_t tmp(sb);
    LowPass(tmp, taps);
    auto nsamps = sb.size();
    for(int j = 0; j < nsamps; ++j) {
        sb[j] = sb[j] - tmp.data[j];
    }
}


void Charge1(sndbfr_48k_t & sb)
{
    sb.set_duration(5);
    auto nsamps = sb.size();
    auto amp = samp_t{0.15};
    for(int j = 0; j < nsamps; ++j) {
        sec_t t = (sec_t)j/sb.FS, t2pi = t*pi2;
        sb[j] = amp*(sin(t2pi*1e3) + sin(t2pi*1.2e3))*sin(t2pi*(t2pi*12));
    }
}


void Bell1(sndbfr_48k_t & sb)
{
    sb.set_duration(2);
    auto nsamps = sb.size();
    auto amp = samp_t{0.15};
    AddDecaySine(sb, 0, 1e3, 0, 0.15, 7);
    AddDecaySine(sb, 0, 0.7e3, 0, 0.15, 11);
    AddDecaySine(sb, 0, 0.3e3, 0, 0.15, 13);
}

void Geiger(sndbfr_48k_t & sb)
{
    auto duration = sec_t{10};
    sb.set_duration(duration);
    for(int j = 0; j < 120; ++j)
        AddDecaySine(sb, rng::frandom(duration), 3e3, 0, 0.15, 500);
}


void Trickle(sndbfr_48k_t & sb)
{
    auto duration = sec_t{10};
    sb.set_duration(duration);
    for(int j = 0; j < 240; ++j)
        AddDecaySine(sb, rng::frandom(duration), rng::frandom(22e3), 0, 0.15, 500);
}

void Crackle(sndbfr_48k_t & sb)
{
    auto duration = sec_t{10};
    sb.set_duration(duration);
    for(int j = 0; j < 240; ++j)
        AddDecaySine(sb, rng::frandom(duration), rng::frandom(2e3, 22e3), 0, 0.15, 5000);
}

void MachineNoise(sndbfr_48k_t & sb)
{
    sb.set_duration(10);
    auto nsamps = sb.size();
    // Machine noise
    auto nwaves = size_t{1024};
    for(auto k = size_t{0}; k < nwaves; ++k)
    {
        freq_t freq = rng::frandom(100, 22e3);
        sec_t phase = rng::signed_random(pi);
//        samp_t amp = 0.5;
        samp_t amp = 10*exp(-freq.f/1e3);
//        samp_t amp = 0.5/(1 + sqr(freq/1e3));
        //AddSine(sb, freq, phase, amp/nwaves);
        for(int j = 0; j < nsamps; ++j) {
            sec_t t = (sec_t)j/sb.FS;
            sb[j] += amp*sin(t*freq.a + phase)/nwaves;
        }
    }
}

void Rushing1(sndbfr_48k_t & sb)
{
    sb.set_duration(15);
    auto nsamps = sb.size();
    auto amp = samp_t{0.5};
    for(int j = 0; j < nsamps; ++j) {
        // auto t = (sec_t)j/sb.FS, t2pi = t*pi2;
        sb[j] = amp*rng::grandom();
    }
//    LowPassMono(sb, sb.FS/0.5e3);
    HighPass(sb, sb.FS/1.5e3);
    for(int j = 0; j < 240; ++j)
        AddDecaySine(sb, rng::frandom(15), rng::frandom(22e3), 0, 0.15, 500);
}


void SineSquareTri(sndbfr_48k_t & sb)
{
    auto amp = samp_t{0.15};
    auto freq = sec_t{1e3};
    auto afreq = freq*pi2;
    sb.set_duration(15);
    auto nsamps = sb.size();
    for(int j = 0; j < nsamps; ++j) {
        auto t = (sec_t)j/sb.FS;
        if(t < 5)
            sb[j] = amp*sin(t*afreq);
        else if(t < 10)
            sb[j] = amp*triangle(t*freq);
        else
            sb[j] = amp*square(t*freq, 0.5);
    }
}


void Warble1(sndbfr_48k_t & sb)
{
    samp_t amp = 0.5;
    sb.set_duration(5);
    auto nsamps = sb.size();
    sec_t t0 = 0, t1 = 0, dt = 1.0/sb.FS;
    for(int j = 0; j < nsamps; ++j) {
        sec_t decay = exp(-t0);
//        freq_t freq = 300 + sin(t0*30*pi2)*100;
        freq_t freq = 300 + sin(t0*30*pi2)*100*decay;
        t0 += dt;
        t1 += (dt*freq.a);
        sb[j] = amp*sin(t1)*decay;
    }
    
//    samp_t amp = 0.5;
//    freq_t freq = 300;
//    sb.set_duration(5);
//    int nsamps = sb.size();
//    for(int j = 0; j < nsamps; ++j) {
//        sec_t t = (sec_t)j/sb.FS;
//        t += sin(t*30*pi2)/30/t;
//        sb[j] = amp*rng::grandom()*sin(t*pi2)*exp(-t));
////        sb[j] = amp*sin(t*freq.a)*exp(-t));
//    }
}


void test_tone(sndbfr_48k_t & sb)
{
    samp_t amp = 0.05;
    sb.set_duration(15);
    auto nsamps = sb.size();
    auto t0 = sec_t{0}, dt = sec_t{1}/sb.FS;
    for(int j = 0; j < nsamps; ++j) {
        t0 += dt;
        auto s = samp_t{0};
        for(int j = 1; j < 11; ++j)
            s += sin(j*333.0*pi2*t0);
        // s = amp*sin(1000.0*pi2*t0);
        // s += amp*sin(2000.0*pi2*t0);
        // s += amp*sin(3000.0*pi2*t0);
        sb[j] = amp*s;
    }
}

void gen_test_signals()
{
    char buf[1024];
    sndbfr_48k_t sb;
    sb.set_duration(60);
    auto nsamps = sb.size();
    
    for(int freq = 1; freq < 22; ++freq)
    {
        for(int j = 0; j < nsamps; ++j) {
            sec_t t0 = (sec_t)j/sb.FS;
            sb[j] = sin(freq*1000.0*pi2*t0)*0.25;
        }
        snprintf(buf, sizeof(buf), "sine%dkHz", freq);
        WriteFLAC(sb, buf);
    }
    
    for(int j = 0; j < nsamps; ++j) {
        sec_t t0 = (sec_t)j/sb.FS;
        sb[j] = 0.0f;
        for(int freq = 1; freq < 22; ++freq)
            sb[j] += sin(freq*1000.0*pi2*t0);
    }
    snprintf(buf, sizeof(buf), "sines1-22kHz");
    sb.normalize_to(0.25f);
    WriteFLAC(sb, buf);
    
    for(int j = 0; j < nsamps; ++j) {
        sec_t t0 = (sec_t)j/sb.FS;
        sb[j] = 0.0f;
        for(int freq = 1; freq < 11; ++freq)
            sb[j] += sin(freq*1000.0*pi2*t0);
    }
    snprintf(buf, sizeof(buf), "sines1-11kHz");
    sb.normalize_to(0.25f);
    WriteFLAC(sb, buf);
    
    for(int j = 0; j < nsamps; ++j)
        sb[j] = rng::signed_random();
    snprintf(buf, sizeof(buf), "flat");
    sb.normalize_to(0.25f);
    WriteFLAC(sb, buf);
    
    for(int j = 0; j < nsamps; ++j)
        sb[j] = rng::grandom();
    snprintf(buf, sizeof(buf), "gaussian");
    sb.normalize_to(0.25f);
    WriteFLAC(sb, buf);
}*/

int main(int argc, const char * argv[])
{
    char buf[1024];
    srandom(time(NULL));
    
    cout << "Generating sound" << endl;
    sec_t startT = GetRealSeconds();
    
    auto && data = std::vector<samp_t>(60*44000);
    auto && sb = sndbfr_48k_t{data};
    
    
    auto && sine_data = std::vector<samp_t>(44000, 0);
    auto && sine_bfr = sndbfr_48k_t{sine_data};
    for(size_t j = 0; j < sine_bfr.size(); ++j)
        sine_bfr[j] = sin(pi2*sine_bfr.t(j));
    
    
    for(size_t j = 0; j < sb.size(); ++j)
        sb[j] = sin(1000.0*pi2*sb.t(j))*0.25;
    snprintf(buf, sizeof(buf), "direct_sine");
    WriteFLAC(sb, buf);
    for(size_t j = 0; j < sb.size(); ++j)
        sb[j] = sine_bfr.mod_lerp(1000.0*sb.t(j))*0.25;
    snprintf(buf, sizeof(buf), "interp_sine");
    WriteFLAC(sb, buf);
    
    
    
    
    // sb.set_duration(5);
    // AddDecaySine(sb, 0, 1000, 0, 0.5, 1);
    // Charge1(sb);
    // Bell1(sb);
    // Geiger(sb);
    // Crackle(sb);
    // MachineNoise(sb);
    // Rushing1(sb);
    // SineSquareTri(sb);
    
    // Warble1(sb);
    // test_tone(sb);
    
//    sb.set_duration(10);
//    size_t nsamps = sb.size();
//    samp_t amp = 6;
//    for(int j = 0; j < nsamps; ++j) {
//        sec_t t = (sec_t)j/sb.FS;
//        sb[j] = amp*rng::grandom();
//    }
//    LowPass(sb, kLeftChannel, sb.FS/(0.05e3));
    
    // gen_test_signals();
    
    cout << "Done in " << (GetRealSeconds() - startT) << " seconds, writing" << endl;
    
    // WriteFLAC(sb, "flacout.flac");
    
    
    
    
    return EXIT_SUCCESS;
}


//******************************************************************************
//    Copyright (c) 2015, Christopher James Huff
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

#ifndef SOUNDBUFFER_H
#define SOUNDBUFFER_H

#include <stdlib.h>
#include <assert.h>
#include <cmath>
#include <cfloat>
#include <fstream>

#include <algorithm>
#include <string>
#include <vector>

#include "FLAC++/metadata.h"
#include "FLAC++/encoder.h"

using samp_t = float;
using sec_t = double;

constexpr double pi = M_PI;
constexpr double pi2 = 2.0*M_PI;

// A frequency value type that auto-computes angular frequency, intended to
// save some multiplications and help keep things consistent.
struct freq_t {
    float f, a;
    freq_t(float fr = 1): f(fr), a(f*pi2) {}
    freq_t & operator=(float fr) {f = fr; a = f*pi2; return *this;}
    freq_t & operator=(const freq_t & fr) {f = fr.f; a = fr.a; return *this;}
};

namespace rng {
// Floating point random: range = [0, 1)
static inline float frandom() {return (float)random()/RAND_MAX;}
// Signed floating point random: range = [-1, 1)
static inline float signed_random() {return 2.0f*(frandom() - 0.5f);}
static inline float signed_random(float scl) {return scl*2.0f*(frandom() - 0.5f);}
// Floating point random: range = [-1, scl)
static inline float frandom(float scl) {return scl*frandom();}
// Floating point random: range = [mn, mx)
static inline float frandom(float mn, float mx) {return (mx - mn)*frandom() + mn;}

// Approximate gaussian random, range = [-1, 1)
static inline float grandom() {
    return (signed_random() + signed_random() + signed_random() +
           signed_random() + signed_random() + signed_random() +
           signed_random() + signed_random() + signed_random() +
           signed_random() + signed_random() + signed_random())/12.0f;
}
} // namespace rng

template<typename T>
auto constexpr sqr(T x) -> T {return x*x;}

template<typename T>
auto constexpr cub(T x) -> T {return x*x*x;}

// Waveforms
static inline samp_t sine(sec_t x) {return sin(x*pi2);}
// Ramp from 0 to 1
static inline samp_t ramp(sec_t x) {return x - floor(x);}
// signed, phase-correct ramp wave
static inline samp_t phramp(sec_t x) {return 2*(x/2 - floor((x - 1.0f)/2.0f)) - 2.0f;}
static inline samp_t triangle(sec_t x) {return -std::min(x - floor(x), ceil(x) - x)*4.0f + 1.0f;}
// 0-1 square wave
static inline samp_t pulse(sec_t x, sec_t duty) {return (ramp(x) < duty)? 1.0f : 0.0f;}
// Signed square wave
static inline samp_t square(sec_t x, sec_t duty) {return (ramp(x) < duty)? -1.0f : 1.0f;}

// template<int _FS>
// struct sndbuf_t {
//     static constexpr int FS = _FS;
    
//     std::vector<samp_t> data;
    
//     sndbuf_t(sec_t len = 0, samp_t val = 0) {
//         set_duration(len);
//         clear(val);
//     }
//     sndbuf_t(const sndbuf_t<_FS> & orig) {
//         data = orig.data;
//     }
    
//     auto set_duration(sec_t secs) {set_size(secs*FS);}
//     auto set_size(size_t nsamps) {data.resize(nsamps);}
    
//     auto size() const -> size_t {return data.size();}
//     auto duration() const -> sec_t {return (sec_t)size()/FS;}
//     auto index(sec_t t) const -> int32_t {return floor(t*FS);}
//     auto t(size_t idx) const -> sec_t {return sec_t(idx)/FS;}
// };

template<int _FS>
struct soundslice_t {
    static constexpr auto FS = _FS;
    using iterator = std::vector<samp_t>::iterator;
    using const_iterator = std::vector<samp_t>::const_iterator;
    
    iterator data_begin;
    iterator data_end;
    
    soundslice_t(std::vector<samp_t> & data):
        data_begin{std::begin(data)},
        data_end{data.end()}
    {}
    soundslice_t(const iterator & b, const iterator & e):
        data_begin{b},
        data_end{e}
    {}
    soundslice_t(const soundslice_t<_FS> & orig) = default;
    
    auto size() const -> size_t {return data_end - data_begin;}
    auto duration() const -> sec_t {return sec_t(size())/FS;}
    
    auto begin() -> iterator {return data_begin;}
    auto end() -> iterator {return data_end;}
    auto begin() const -> const_iterator {return data_begin;}
    auto end() const -> const_iterator {return data_end;}
    
    auto slice(size_t start_idx, size_t end_idx) -> soundslice_t<FS> {
        return {begin() + start_idx, begin() + end_idx};
    }
    auto slice(size_t start_idx) -> soundslice_t<FS> {
        return {begin() + start_idx, end()};
    }
    auto slice(sec_t tstart, sec_t tend) -> soundslice_t<FS> {
        return {begin() + index(tstart), begin() + index(tend)};
    }
    auto slice(sec_t tstart) -> soundslice_t<FS> {
        return {begin() + index(tstart), end()};
    }
    
    auto index(sec_t t) const -> int32_t {return floor(t*FS);}
    auto t(size_t idx) const -> sec_t {return sec_t(idx)/FS;}
    
    auto operator[](size_t idx) -> samp_t & {return *(data_begin + idx);}
    auto operator[](size_t idx) const -> const samp_t & {return *(data_begin + idx);}
    
    auto operator()(sec_t t) -> samp_t & {return *(data_begin + index(t));}
    auto operator()(sec_t t) const -> const samp_t & {return *(data_begin + index(t));}
    
    auto lerp(sec_t t) -> samp_t {
        auto integ = sec_t{};
        auto frac = modf(t*FS, &integ);
        auto idx = size_t(integ);
        return (*this)[idx]*(sec_t{1} - frac) + (*this)[idx + 1]*frac;
    }
    auto bounded_lerp(sec_t t) -> samp_t {
        auto integ = sec_t{};
        auto frac = modf(t*FS, &integ);
        auto idx = size_t(integ);
        if(idx < size())
            return (*this)[idx]*(sec_t{1} - frac) + (*this)[idx + 1]*frac;
        else
            return samp_t{0};
    }
    auto mod_lerp(sec_t t) -> samp_t {
        auto integ = sec_t{};
        auto frac = modf(t*FS, &integ);
        auto idx = size_t(integ)%size();
        return (*this)[idx]*(sec_t{1} - frac) + (*this)[idx + 1]*frac;
    }
    
    void copy(const soundslice_t<_FS> & rhs, size_t start_idx = 0) {
        auto && k = begin(rhs);
        for(auto && s: slice(start_idx, std::min(size(), rhs.size())))
            s = *k++;
    }
    void add(const soundslice_t<_FS> & rhs, size_t start_idx = 0) {
        auto && k = begin(rhs);
        for(auto && s: slice(start_idx, std::min(size(), rhs.size())))
            s += *k++;
    }
    void mul(const soundslice_t<_FS> & rhs, size_t start_idx = 0) {
        auto && k = begin(rhs);
        for(auto && s: slice(start_idx, std::min(size(), rhs.size())))
            s *= *k++;
    }
    
    
    auto clear(samp_t val, sec_t tstart, sec_t tend) {
        slice(tstart, tend).clear();
    }
    auto clear(samp_t val) {
        for(auto && s: *this)
            s = val;
    }
    
    auto add(samp_t val) {
        for(auto && s: *this)
            s += val;
    }
    
    auto mul(samp_t val) {
        for(auto && s: *this)
            s *= val;
    }
    
    
    auto calc_min_max(samp_t & mn, samp_t & mx) const {
        mn = FLT_MAX;
        mx = FLT_MIN;
        for(auto && s: *this)
        {
            if(s < mn)
                mn = s;
            if(s > mx)
                mx = s;
        }
    }
    
    auto calc_max_amp(samp_t & mx) const {
        mx = 0;
        for(auto && s: *this)
        {
            if(-s > mx)
                mx = -s;
            if(s > mx)
                mx = s;
        }
    }
    
    auto normalize_to(samp_t max_amp) {
        samp_t mx;
        calc_max_amp(mx);
        samp_t scl = max_amp/mx;
        for(auto && s: *this)
            s *= scl;
    }
    
    auto debug() const
    {
        samp_t mn, mx;
        calc_min_max(mn, mx);
        std::cout << "duration: " << duration() << ", size: " << size() << std::endl;
        std::cout << "min sample: " << mn << ", max: " << mx << std::endl;
    }
};


// enum ChannelIdx {
//     kLeftChannel = 0,
//     kFrontLeftChannel = kLeftChannel,
//     kRightChannel = 1,
//     kFrontRightChannel = kRightChannel,
//     kFrontCenterChannel = 2,
//     kLFE_Channel = 3,
//     kBackLeftChannel = 4,
//     kBackRightChannel = 5,
//     kMaxChannels
// };

template<typename bfr_t>
void write_s32(bfr_t & buf, const std::string & path) {
    buf.debug();
    char fname[1024];
    snprintf(fname, sizeof(fname), "%s.s32", path.c_str());
    auto && out = std::ofstream(fname, std::ofstream::binary);
    out.write((char *)&buf[0], buf.size()*sizeof(float));
    
    char cmd[1024];
    snprintf(cmd, sizeof(cmd), "sox -r 48k -c 1 %s -e float %s.wav",
        fname, path.c_str());
    printf("%s\n", cmd);
    system(cmd);
}

template<typename bfr_t>
void WriteFLAC(bfr_t & buf, const std::string & path) {
    buf.debug();
    FLAC::Encoder::File enc;
    
    auto channels = 1;
    bool ok = true;
    ok &= enc.set_verify(true);
    ok &= enc.set_compression_level(5);
    ok &= enc.set_channels(channels);
    ok &= enc.set_bits_per_sample(16);
    ok &= enc.set_sample_rate(bfr_t::FS);
    ok &= enc.set_total_samples_estimate(buf.size()*channels);
    
    if(!ok) {
        std::cerr << "Error: couldn't configure encoder" << std::endl;
        return;
    }
    
    auto && fname = path + ".flac";
    FLAC__StreamEncoderInitStatus init_status = enc.init(fname.c_str());
    if(init_status != FLAC__STREAM_ENCODER_INIT_STATUS_OK) {
        std::cerr << "Error: couldn't init encoder" << std::endl;
        return;
    }
    
    std::vector<FLAC__int32> data(channels*buf.size());
    FLAC__int32 * pcmData[8];
    for(size_t ch = 0; ch < channels; ++ch)
    {
        pcmData[ch] = &(data[ch*buf.size()]);
        FLAC__int32 * dst = pcmData[ch];
        for(auto && samp: buf)
            *dst++ = (FLAC__int32)(samp*0x7FFF);
            // *dst++ = (FLAC__int32)(samp*0x007FFFFF);
    }
    
    if(!enc.process(pcmData, buf.size())) {
        std::cerr << "Error while encoding" << std::endl;
    }
    
    enc.finish();
}

typedef soundslice_t<48000> sndbfr_48k_t;

typedef soundslice_t<192000> sndbfr_192k_t;

#endif // SOUNDBUFFER_H

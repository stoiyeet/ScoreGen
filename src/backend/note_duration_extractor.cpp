#define _USE_MATH_DEFINES
#include "note_duration_extractor.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <iostream>

// Constants for pitch detection (avoids magic numbers)
const double MIN_PITCH_HZ = 50.0;    // Lowest detectable pitch (e.g., for bass notes)
const double MAX_PITCH_HZ = 2000.0;  // Highest detectable pitch
const double PITCH_CONFIDENCE_THRESHOLD = 0.8;  // Autocorrelation must be at least this strong
const double RMS_THRESHOLD = 0.001;  // Minimum RMS for pitch detection

//
// Function: detectPitch
// ---------------------
// Uses autocorrelation to estimate the dominant pitch (in Hz) from a windowed frame.
// Returns 0 if no clear pitch is detected.
//
double detectPitch(const std::vector<double>& frame, int sampleRate) {
    size_t N = frame.size();
    if (N == 0) return 0.0;
    
    // Compute autocorrelation
    std::vector<double> autocorr(N, 0.0);
    double maxAutocorr = 0.0;
    for (size_t lag = 0; lag < N; ++lag) {
        for (size_t i = 0; i < N - lag; ++i) {
            autocorr[lag] += frame[i] * frame[i + lag];
        }
        if (lag == 0) {
            maxAutocorr = autocorr[0];
        }
    }
    
    // Skip zero lag and find the lag with maximum autocorrelation
    size_t maxLag = 0;
    double maxVal = 0.0;
    size_t minLag = static_cast<size_t>(std::ceil(sampleRate / MAX_PITCH_HZ));
    size_t maxLagLimit = static_cast<size_t>(std::floor(sampleRate / MIN_PITCH_HZ));
    
    for (size_t lag = minLag; lag <= maxLagLimit && lag < N; ++lag) {
        if (autocorr[lag] > maxVal) {
            maxVal = autocorr[lag];
            maxLag = lag;
        }
    }
    
    // Check if the autocorrelation is strong enough relative to zero lag
    if (maxLag > 0 && (maxVal / maxAutocorr) > PITCH_CONFIDENCE_THRESHOLD) {
        return static_cast<double>(sampleRate) / maxLag;
    }
    return 0.0;
}

//
// Function: frequencyToNoteString
// -------------------------------
// Converts a frequency (Hz) to a musical note name (e.g., "A4", "C#3").
// Uses A4 = 440 Hz as reference. Frequencies <= 0 return "Rest".
//
std::string frequencyToNoteString(double frequency) {
    if (frequency <= 0.0) return "Rest";
    
    // Constants for note calculation
    const double A4_FREQUENCY = 440.0;
    const int A4_OCTAVE = 4;
    const std::vector<std::string> NOTE_NAMES = {
        "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"
    };
    
    // Calculate number of semitones from A4
    double semitonesFromA4 = 12.0 * std::log2(frequency / A4_FREQUENCY);
    int roundedSemitones = static_cast<int>(std::round(semitonesFromA4));
    
    // Calculate octave and note index
    int noteIndex = (roundedSemitones + 9) % 12;  // +9 because A is 9th in NOTE_NAMES
    int octaveOffset = (roundedSemitones + 9) / 12;
    int octave = A4_OCTAVE + octaveOffset;
    
    // Adjust for negative indices
    if (noteIndex < 0) noteIndex += 12;
    
    return NOTE_NAMES[noteIndex] + std::to_string(octave);
}


// Internal structure to hold note segment information.
struct NoteSegment {
    std::string note;
    int startFrame; // First frame index of the note.
    int endFrame;   // Last frame index of the note.
};

std::string determineNoteType(float noteDuration, int bpm) {
    // Predefined note durations
    std::map<std::string, float> note_durations = {
        {"sixteenth", 0.25},
        {"eighth", 0.5},
        {"quarter", 1.0},
        {"dotted quarter", 1.5},
        {"half", 2.0},
        {"dotted half", 3.0},
        {"whole", 4.0}
    };

    float beatDuration = 60.0f / bpm;
    float beatsPerNote = noteDuration / beatDuration;

    auto closest = min_element(
        note_durations.begin(),
        note_durations.end(),
        [beatsPerNote](const auto& a, const auto& b) {
            return abs(a.second - beatsPerNote) < abs(b.second - beatsPerNote);
        }
    );
    return closest->first;
}

//
// Function: extract_note_durations
// --------------------------------
// Processes the WAV file and returns a vector of Note objects that contain
// start time, end time, pitch (as a note string), and note type (set to "unknown").
// Now also performs a simple onset detection: if an onset is detected in the middle
// of a note segment, that segment is split into multiple notes.
//
std::vector<Note> extract_note_durations(const char* infilename, int bpm) {
    std::vector<Note> notes;
    std::vector<double> audio;
    int sampleRate;
    if (!readWav(infilename, audio, sampleRate)) {
        std::cerr << "Error reading WAV file.\n";
        return notes;
    }
    
    // Analysis parameters for pitch detection.
    const int FRAME_SIZE = 2048;  // Window size for pitch detection
    const int HOP_SIZE = 512;     // Hop size for pitch detection
    
    int totalSamples = static_cast<int>(audio.size());
    int numFrames = (totalSamples >= FRAME_SIZE) ? ((totalSamples - FRAME_SIZE) / HOP_SIZE + 1) : 0;
    
    std::vector<double> window = hanningFunction(FRAME_SIZE);
    std::vector<double> pitchEstimates(numFrames, 0.0);
    std::vector<double> frameRMS(numFrames, 0.0); // RMS energy per pitch frame
    
    // Compute overall RMS to make thresholds relative
    double totalSumSq = 0.0;
    for (double sample : audio) {
        totalSumSq += sample * sample;
    }
    double overallRMS = std::sqrt(totalSumSq / audio.size());
    double relativeRMSThreshold = overallRMS * 0.1;  // 10% of overall RMS
    
    // Compute pitch estimates and RMS using the larger window.
    for (int frame = 0; frame < numFrames; frame++) {
        int start = frame * HOP_SIZE;
        std::vector<double> frameBuffer(FRAME_SIZE);
        double sumSq = 0.0;
        for (int n = 0; n < FRAME_SIZE; n++) {
            frameBuffer[n] = audio[start + n] * window[n];
            sumSq += frameBuffer[n] * frameBuffer[n];
        }
        double rms = std::sqrt(sumSq / FRAME_SIZE);
        frameRMS[frame] = rms;
        if (rms < std::max(RMS_THRESHOLD, relativeRMSThreshold)) {
            pitchEstimates[frame] = 0.0;
        } else {
            double pitch = detectPitch(frameBuffer, sampleRate);
            pitchEstimates[frame] = pitch;
        }
    }
        
    // Segment frames into note and rest segments.
    double tolerance = 0.05;         // allow ~4% pitch variation within a note
    double minNoteDuration = 60.0 / (bpm * 4); // minimum segment duration in seconds
    bool inSegment = false;
    bool isNoteSegment = false;      // true if current segment is a note, false if a rest
    double currentPitch = 0.0;       // used if in a note segment
    int segmentStartFrame = 0;
    std::vector<NoteSegment> segments;
    
    for (int i = 0; i < numFrames; i++) {
        double pitch = pitchEstimates[i];
        bool isRest = (pitch == 0);
        if (!inSegment) {
            inSegment = true;
            segmentStartFrame = i;
            isNoteSegment = !isRest;
            if (isNoteSegment)
                currentPitch = pitch;
        } else {
            if (isNoteSegment) {
                // End the note segment if a rest occurs or pitch deviates too much.
                if (isRest || std::abs(pitch - currentPitch) / currentPitch > tolerance) {
                    int segmentEndFrame = i - 1;
                    double startTime = segmentStartFrame * HOP_SIZE / static_cast<double>(sampleRate);
                    double endTime = (segmentEndFrame * HOP_SIZE + FRAME_SIZE) / static_cast<double>(sampleRate);
                    double duration = endTime - startTime;
                    if (duration >= minNoteDuration) {
                        segments.push_back({frequencyToNoteString(currentPitch), segmentStartFrame, segmentEndFrame});
                    }
                    // Start a new segment.
                    inSegment = true;
                    segmentStartFrame = i;
                    isNoteSegment = !isRest;
                    if (isNoteSegment)
                        currentPitch = pitch;
                }
            } else {
                // In a rest segment, if a note is encountered, end the rest segment.
                if (!isRest) {
                    int segmentEndFrame = i - 1;
                    double startTime = segmentStartFrame * HOP_SIZE / static_cast<double>(sampleRate);
                    double endTime = (segmentEndFrame * HOP_SIZE + FRAME_SIZE) / static_cast<double>(sampleRate);
                    double duration = endTime - startTime;
                    if (duration >= minNoteDuration) {
                        segments.push_back({"Rest", segmentStartFrame, segmentEndFrame});
                    }
                    inSegment = true;
                    segmentStartFrame = i;
                    isNoteSegment = true;
                    currentPitch = pitch;
                }
            }
        }
    }
    if (inSegment) {
        int segmentEndFrame = numFrames - 1;
        double startTime = segmentStartFrame * HOP_SIZE / static_cast<double>(sampleRate);
        double endTime = (segmentEndFrame * HOP_SIZE + FRAME_SIZE) / static_cast<double>(sampleRate);
        double duration = endTime - startTime;
        if (duration >= minNoteDuration) {
            if (isNoteSegment) {
                segments.push_back({frequencyToNoteString(currentPitch), segmentStartFrame, segmentEndFrame});
            } else {
                segments.push_back({"Rest", segmentStartFrame, segmentEndFrame});
            }
        }
    }

    // Merge adjacent segments with the same note if the gap is small.
    double mergeThreshold = 0.05; // seconds
    std::vector<NoteSegment> mergedSegments;
    if (!segments.empty()) {
        mergedSegments.push_back(segments[0]);
        for (size_t i = 1; i < segments.size(); i++) {
            NoteSegment &prev = mergedSegments.back();
            NoteSegment &curr = segments[i];
            double prevEndTime = (prev.endFrame * HOP_SIZE + FRAME_SIZE) / static_cast<double>(sampleRate);
            double currStartTime = (curr.startFrame * HOP_SIZE) / static_cast<double>(sampleRate);
            if (prev.note == curr.note && (currStartTime - prevEndTime) < mergeThreshold) {
                prev.endFrame = curr.endFrame;
            } else {
                mergedSegments.push_back(curr);
            }
        }
    }
    
    // --- Onset Detection with a Smaller Window ---
    const int ONSET_FRAME_SIZE = 512;  // Smaller window for onset detection
    const int ONSET_HOP_SIZE = 256;    // Higher time resolution
    int numOnsetFrames = (totalSamples >= ONSET_FRAME_SIZE) ?
                         ((totalSamples - ONSET_FRAME_SIZE) / ONSET_HOP_SIZE + 1) : 0;
    std::vector<double> onsetRMS(numOnsetFrames, 0.0);
    for (int i = 0; i < numOnsetFrames; i++) {
        int start = i * ONSET_HOP_SIZE;
        double sumSq = 0.0;
        for (int n = 0; n < ONSET_FRAME_SIZE; n++) {
            sumSq += audio[start + n] * audio[start + n];
        }
        onsetRMS[i] = std::sqrt(sumSq / ONSET_FRAME_SIZE);
    }
    
    // Collect onset times - consider first frame as potential onset if above threshold
    std::vector<double> onsetTimes;
    double onsetThreshold = overallRMS * 0.2;  // Relative threshold
    if (onsetRMS[0] > onsetThreshold) {
        onsetTimes.push_back(0.0);  // Consider start as onset if significant energy
    }
    for (int i = 1; i < numOnsetFrames; i++) {
        double diff = onsetRMS[i] - onsetRMS[i - 1];
        if (diff > onsetThreshold && onsetRMS[i] > relativeRMSThreshold) {
            double T = (i * ONSET_HOP_SIZE) / static_cast<double>(sampleRate);
            onsetTimes.push_back(T);
        }
    }
    
    // --- Split Note Segments at Onsets ---
    // For each merged note segment (ignoring rests), check if any onset (from the small-window analysis)
    // occurs within its time boundaries. If so, split the segment at the corresponding pitch-frame indices.
    std::vector<NoteSegment> finalSegments;
    for (const auto &seg : mergedSegments) {
        if (seg.note != "Rest") {
            double segStartTime = seg.startFrame * HOP_SIZE / static_cast<double>(sampleRate);
            double segEndTime = (seg.endFrame * HOP_SIZE + FRAME_SIZE) / static_cast<double>(sampleRate);
            std::vector<int> splitPoints;
            // For each detected onset time, if it falls inside the segment, map it to a pitch-frame index.
            for (double T : onsetTimes) {
                if (T > segStartTime && T < segEndTime) {
                    int pitchFrameIndex = static_cast<int>(std::round(T * sampleRate / HOP_SIZE));
                    if (pitchFrameIndex > seg.startFrame && pitchFrameIndex <= seg.endFrame)
                        splitPoints.push_back(pitchFrameIndex);
                }
            }
            if (splitPoints.empty()) {
                finalSegments.push_back(seg);
            } else {
                int currentStart = seg.startFrame;
                for (int pf : splitPoints) {
                    int segmentEnd = pf - 1;
                    double startTime = currentStart * HOP_SIZE / static_cast<double>(sampleRate);
                    double endTime = (segmentEnd * HOP_SIZE + FRAME_SIZE) / static_cast<double>(sampleRate);
                    double duration = endTime - startTime;
                    if (duration >= minNoteDuration)
                        finalSegments.push_back({seg.note, currentStart, segmentEnd});
                    currentStart = pf;
                }
                double startTimeFinal = currentStart * HOP_SIZE / static_cast<double>(sampleRate);
                double endTimeFinal = (seg.endFrame * HOP_SIZE + FRAME_SIZE) / static_cast<double>(sampleRate);
                double durationFinal = endTimeFinal - startTimeFinal;
                if (durationFinal >= minNoteDuration)
                    finalSegments.push_back({seg.note, currentStart, seg.endFrame});
            }
        } else {
            finalSegments.push_back(seg);
        }
    }
    
    // Convert final segments to Note objects.
    for (const auto &seg : finalSegments) {
        double startTime = seg.startFrame * HOP_SIZE / static_cast<double>(sampleRate);
        double endTime = (seg.endFrame * HOP_SIZE + FRAME_SIZE) / static_cast<double>(sampleRate);
        std::string noteType = determineNoteType((endTime - startTime), bpm);
        notes.push_back({static_cast<float>(startTime), static_cast<float>(endTime), seg.note, noteType});
        std::cout << "Note: " << seg.note << " | Start Time: " << startTime 
                  << " s | End Time: " << endTime << " s | Type: " << noteType << "\n";
    }
    return notes;
}

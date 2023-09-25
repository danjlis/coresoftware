#ifndef FUN4ALLRAW_SINGLETPCINPUT_H
#define FUN4ALLRAW_SINGLETPCINPUT_H

#include "SingleStreamingInput.h"

#include "InttPool.h"

#include <array>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

class Eventiterator;
class Fun4AllEvtInputPoolManager;
class InttRawHit;
class Packet;

class SingleTPCInput : public SingleStreamingInput
{
 public:
  explicit SingleTPCInput(const std::string &name);
  ~SingleTPCInput() override;
  void FillPool(const unsigned int nevents = 1) override;
  void CleanupUsedPackets(const uint64_t bclk) override;
  bool CheckPoolDepth(const uint64_t bclk) override;
  void ClearCurrentEvent() override;
  void Print(const std::string &what = "ALL") const override;

 private:
  Packet **plist = nullptr;
  InttPool pool;
  unsigned int m_NumSpecialEvents = 0;
  std::array<uint64_t, 14> m_PreviousClock{};
  std::array<uint64_t, 14> m_Rollover{};
  std::map<uint64_t, std::set<int>> m_BeamClockFEE;
  std::map<uint64_t, std::vector<InttRawHit *>> m_InttRawHitMap;
  std::map<int, uint64_t> m_FEEBclkMap;
  std::set<uint64_t> m_BclkStack;
};

#endif

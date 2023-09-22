#ifndef FUN4ALLRAW_SINGLESTREAMINGINPUT_H
#define FUN4ALLRAW_SINGLESTREAMINGINPUT_H

#include <fun4all/Fun4AllBase.h>
#include <fun4all/InputFileHandler.h>

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

class SingleStreamingInput : public Fun4AllBase, public InputFileHandler
{
 public:
  explicit SingleStreamingInput(const std::string &name, Fun4AllEvtInputPoolManager *inman);
  explicit SingleStreamingInput(const std::string &name);
  ~SingleStreamingInput() override;
  virtual Eventiterator *GetEventIterator() { return m_EventIterator; }
  virtual void FillPool(const unsigned int nevents = 1);
  virtual void RunNumber(const int runno) { m_RunNumber = runno; }
  virtual int RunNumber() const { return m_RunNumber; }
  virtual int fileopen(const std::string &filename) override;
  virtual int fileclose() override;
  virtual int AllDone() const { return m_AllDone; }
  virtual void AllDone(const int i) { m_AllDone = i; }
  virtual void EventNumberOffset(const int i) { m_EventNumberOffset = i; }
  virtual void Print(const std::string &what = "ALL") const override;
  virtual void CleanupUsedPackets(const uint64_t bclk);
  virtual bool CheckPoolDepth(const uint64_t bclk);
  virtual void ClearCurrentEvent();
  virtual Eventiterator *GetEventiterator() const {return m_EventIterator;}
  virtual Fun4AllEvtInputPoolManager *InputManager() {return m_InputMgr;}
  virtual void InputManager(Fun4AllEvtInputPoolManager *in) {m_InputMgr = in;}

 private:
  Eventiterator *m_EventIterator = nullptr;
  Fun4AllEvtInputPoolManager *m_InputMgr = nullptr;
  Packet **plist = nullptr;
  unsigned int m_NumSpecialEvents = 0;
  unsigned int m_EventNumberOffset = 1;  // packet event counters start at 0 but we start with event number 1
  int m_RunNumber = 0;
  int m_EventsThisFile = 0;
  int m_AllDone = 0;
  std::array<uint64_t, 14> m_PreviousClock{};
  std::array<uint64_t, 14> m_Rollover{};
  std::map<uint64_t, std::set<int>> m_BeamClockFEE;
  std::map<uint64_t, std::vector<InttRawHit *>> m_InttRawHitMap;
  std::map<int, uint64_t> m_FEEBclkMap;
  std::set<uint64_t> m_BclkStack;
};

#endif

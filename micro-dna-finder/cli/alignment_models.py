from dataclasses import dataclass

@dataclass(slots=True)
class AlignmentDetails:
    seq: str
    sc_len: int = 0
    count: int = 1
    
    def update_count(self):
        self.count += 1

@dataclass(slots=True)
class AlignerDetails:
    alignment: str
    score: float


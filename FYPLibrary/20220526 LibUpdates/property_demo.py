class Student:
    def __init__(self):
        self._score = 0

    @property
    def score(self):
        return self._score

    @score.setter
    def score(self, s):
        if 0 <= s <= 100:
            self._score = s
        else:
            raise ValueError('The score must be between 0 ~ 100!')

    @score.deleter
    def score(self):
        del self._score

if __name__ == "__main__":
    James = Student()
    James.score = 45
    print(f"James scored {James.score}.")
    Yang = Student()
    Yang.score = 999
    # ValueError: The score must be between 0 ~ 100!
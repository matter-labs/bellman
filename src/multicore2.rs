use std::collections::VecDeque;
use std::sync::Arc;

use parking_lot::{Condvar, Mutex};

use crossbeam::channel::{self, Receiver};

type Task = Box<dyn FnOnce() + Send + 'static>;

#[derive(Default)]
struct TaskQueue(VecDeque<Task>, VecDeque<Task>);

#[derive(Default)]
struct State(Mutex<TaskQueue>, Condvar);

#[derive(Clone)]
pub struct Worker {
    pub(crate) cpus: usize,
    shared_state: Arc<State>,
}

impl Worker {
    fn spawn_in_pool(shared_state: Arc<State>) {
        std::thread::spawn(move || {
            let State(ref mutex, ref cv) = &*shared_state;
            let mut guard = mutex.lock();

            loop {
                while let Some(task) = guard.0.pop_front() {
                    drop(guard);
                    task();
                    guard = mutex.lock();
                }

                match guard.1.pop_front() {
                    Some(task) => {
                        drop(guard);
                        task();
                        guard = mutex.lock();
                    }
                    None => cv.wait(&mut guard),
                }
            }
        });
    }

    fn start_workers(cpus: usize, shared_state: &Arc<State>) {
        for _ in 0..cpus {
            Self::spawn_in_pool(shared_state.clone());
        }
    }

    pub fn new_with_cpus(cpus: usize) -> Worker {
        assert!(cpus > 0);

        let shared_state = Arc::default();
        Self::start_workers(cpus, &shared_state);
        Worker { cpus, shared_state }
    }

    pub fn new() -> Worker {
        Self::new_with_cpus(num_cpus::get())
    }

    pub fn compute<F, T, E>(&self, f: F) -> Receiver<Result<T, E>>
    where
        F: FnOnce() -> Result<T, E> + Send + 'static,
        T: Send + 'static,
        E: Send + 'static,
    {
        let State(ref mutex, ref cv) = &*self.shared_state;
        let (sender, receiver) = channel::bounded(1);

        let boxed_fn = Box::new(move || {
            let result = f();
            sender.send(result).unwrap();
        });

        {
            let mut guard = mutex.lock();
            guard.0.push_back(boxed_fn);
        }
        cv.notify_one();

        receiver
    }

    pub fn compute_background<F, T, E>(&self, f: F) -> Receiver<Result<T, E>>
    where
        F: FnOnce() -> Result<T, E> + Send + 'static,
        T: Send + 'static,
        E: Send + 'static,
    {
        let State(ref mutex, ref cv) = &*self.shared_state;
        let (sender, receiver) = channel::bounded(1);

        let boxed_fn = Box::new(move || {
            let result = f();
            sender.send(result).unwrap();
        });

        {
            let mut guard = mutex.lock();
            guard.1.push_back(boxed_fn);
        }
        cv.notify_one();

        receiver
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_compute() {
        let worker = Worker::new();

        let receiver: Receiver<Result<u64, ()>> = worker.compute(|| {
            let mut i = 0;
            for j in 1..1000000 {
                i += j;
            }
            Ok(i)
        });

        let result = receiver.recv().unwrap();

        assert_eq!(result, Ok(499999500000));
    }
}
